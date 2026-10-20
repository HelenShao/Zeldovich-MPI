/*
 * Kernel_MPI: MPI/C port of 1rank-run-cosm500/kernel_MPI.py
 *
 * Reads the Gaussian potential y-slices dumped by Zeldovich_MPI (DUMP_GAUSSIAN_D),
 *     <BIN_DIR>/Phi_gaussian_rank{r}_y{y}.bin
 * computes the f_NL kernel selected in the par2 file, and writes
 *     <OUT_DIR>/phi_output_rank{r}_y{y}.bin    (default OUT_DIR = <BIN_DIR>/bin_files_output)
 * for the LOAD_D_FROM_FILE run. Must be launched with the same number of ranks and
 * NumZRanks as the Zeldovich_MPI run that wrote the dumps (rank r reads/writes the
 * y_primary slices it generated).
 *
 * Usage:  mpiexec -n R Kernel_MPI [param.par2] [--bin-dir DIR] [--out-dir DIR]
 *
 * Par2 keys:  BoxSize, NumZRanks, ZD_f_NL_local, ZD_f_NL_equil, ZD_f_NL_orthoCMB,
 *             ZD_f_NL_orthoLSS, optional KernelPadM (padded grid size M)
 *
 * Environment:
 *   PARAM_FILE, BIN_DIR, OUT_DIR    same meaning as in kernel_MPI.py (CLI args override)
 *   KERNEL_PYTHON_COMPAT=1 (default) default M = 2N, as in kernel_MPI.py
 *                       =0           default M = smallest valid M >= 3N/2 (alias-free, smaller)
 *   KERNEL_PAD_M=<M>                 explicit padded size (overrides KernelPadM and the default)
 *   KERNEL_ZERO_NYQUIST=1            zero the kx/ky/kz = N/2 modes of phi before the kernels
 *   KERNEL_ONLY_TERM=<A1|C1|...>     debug: run only this shape with coef 1 (output = J*Q)
 *   KERNEL_FORCE_STAGED_A2A=1        test: always use the INT_MAX-safe staged Alltoallv path
 *
 * Stages per kernel term: see kernel_pass.h.
 */

#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <sys/stat.h>

#include "config.h"
#include "precision.h"
#include "types.h"
#include "mpi_topology.h"

#include "kernel_grid.h"
#include "kernel_io.h"
#include "kernel_pass.h"
#include "kernel_terms.h"

#ifndef USE_DOUBLE_PRECISION
#error "Kernel_MPI requires double precision: configure meson with -Ddouble_precision=true"
#endif

// ====================================================================================
// par2 parsing (port of parse_par2 in kernel_MPI.py, numeric values only)
// ====================================================================================

typedef struct {
    double box_size;
    double num_z_ranks;
    double fnl[4];       // local, equil, orthoCMB, orthoLSS
    double pad_M;        // KernelPadM, or -1
    int found[7];        // BoxSize, NumZRanks, 4 x f_NL, KernelPadM
} KernelParams;

static const char *param_keys[7] = {
    "BoxSize", "NumZRanks", "ZD_f_NL_local", "ZD_f_NL_equil",
    "ZD_f_NL_orthoCMB", "ZD_f_NL_orthoLSS", "KernelPadM"
};

static char *strip(char *s)
{
    while (isspace((unsigned char)*s)) s++;
    char *e = s + strlen(s);
    while (e > s && isspace((unsigned char)e[-1])) *--e = '\0';
    return s;
}

static int parse_par2(const char *path, KernelParams *kp)
{
    FILE *fp = fopen(path, "r");
    if (fp == NULL) return 1;
    memset(kp, 0, sizeof(*kp));
    kp->pad_M = -1.0;

    char line[4096];
    while (fgets(line, sizeof(line), fp) != NULL) {
        char *hash = strchr(line, '#');
        if (hash) *hash = '\0';
        char *eq = strchr(line, '=');
        if (eq == NULL) continue;
        *eq = '\0';
        char *key = strip(line);
        char *val = strip(eq + 1);
        size_t vlen = strlen(val);
        if (vlen >= 2 && (val[0] == '"' || val[0] == '\'') && val[vlen - 1] == val[0]) {
            val[vlen - 1] = '\0';
            val = strip(val + 1);
        }

        for (int i = 0; i < 7; i++) {
            if (strcmp(key, param_keys[i]) != 0) continue;
            char *end = NULL;
            errno = 0;
            double v = strtod(val, &end);
            if (end == val || errno != 0) break;   // non-numeric: ignore, like the Python fallback
            double *dst = (i == 0) ? &kp->box_size
                        : (i == 1) ? &kp->num_z_ranks
                        : (i == 6) ? &kp->pad_M
                        : &kp->fnl[i - 2];
            *dst = v;
            kp->found[i] = 1;
        }
    }
    fclose(fp);
    return 0;
}

static int env_int(const char *name, int dflt)
{
    const char *v = getenv(name);
    return (v != NULL && v[0] != '\0') ? atoi(v) : dflt;
}

static double peak_rss_gb(void)
{
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    return ru.ru_maxrss / 1.0e6;   // ru_maxrss is in kB on Linux
}

#define FAIL(...) do { if (world_rank == 0) fprintf(stderr, "[Kernel_MPI] ERROR: " __VA_ARGS__); \
                       MPI_Abort(MPI_COMM_WORLD, 1); } while (0)

int main(int argc, char **argv)
{
    int provided = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    int world_rank = 0, num_ranks = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);
    double t_start = MPI_Wtime();

    // ==================================================================================
    // Arguments / environment
    // ==================================================================================
    const char *param_file = getenv("PARAM_FILE") ? getenv("PARAM_FILE") : ".runtime/param.par2";
    const char *bin_dir = getenv("BIN_DIR") ? getenv("BIN_DIR") : "bin_files";
    const char *out_dir_env = getenv("OUT_DIR");
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--bin-dir") == 0 && i + 1 < argc) {
            bin_dir = argv[++i];
        } else if (strcmp(argv[i], "--out-dir") == 0 && i + 1 < argc) {
            out_dir_env = argv[++i];
        } else if (argv[i][0] != '-') {
            param_file = argv[i];
        } else {
            FAIL("unknown argument %s\n", argv[i]);
        }
    }
    char out_dir[4096];
    if (out_dir_env != NULL && out_dir_env[0] != '\0') {
        snprintf(out_dir, sizeof(out_dir), "%s", out_dir_env);
    } else {
        snprintf(out_dir, sizeof(out_dir), "%s/bin_files_output", bin_dir);
    }

    const int compat = env_int("KERNEL_PYTHON_COMPAT", 1);
    const int zero_nyquist = env_int("KERNEL_ZERO_NYQUIST", 0);
    const int pad_M_env = env_int("KERNEL_PAD_M", -1);
    const char *only_term = getenv("KERNEL_ONLY_TERM");

    // ==================================================================================
    // Rank 0: parameters, N, dumped rank count -> broadcast
    // ==================================================================================
    KernelParams kp;
    memset(&kp, 0, sizeof(kp));
    int hdr[3] = {0, -1, -1};   // parse status, N, max dumped rank
    if (world_rank == 0) {
        hdr[0] = parse_par2(param_file, &kp);
        char path[4096];
        snprintf(path, sizeof(path), "%s/Phi_gaussian_rank0_y0.bin", bin_dir);
        hdr[1] = kernel_detect_N(path);
        hdr[2] = kernel_max_dumped_rank(bin_dir);
    }
    MPI_Bcast(hdr, 3, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&kp, sizeof(kp), MPI_BYTE, 0, MPI_COMM_WORLD);

    if (hdr[0] != 0) FAIL("cannot read parameter file %s\n", param_file);
    if (!kp.found[0] || !kp.found[1]) FAIL("%s must define BoxSize and NumZRanks\n", param_file);
    const int N = hdr[1];
    if (N <= 0) FAIL("%s/Phi_gaussian_rank0_y0.bin missing or not N*N complex doubles\n", bin_dir);
    if (N % 2 != 0) FAIL("N=%d must be even\n", N);
    if (hdr[2] + 1 != num_ranks) {
        FAIL("%s holds dumps from %d ranks but Kernel_MPI runs on %d; launch with the same rank count "
             "as the Zeldovich_MPI run (or clear old dumps)\n", bin_dir, hdr[2] + 1, num_ranks);
    }
    const int num_z_ranks = (int)kp.num_z_ranks;
    if (num_z_ranks <= 0 || num_ranks % num_z_ranks != 0) {
        FAIL("num_ranks=%d is not divisible by NumZRanks=%d\n", num_ranks, num_z_ranks);
    }
    if (num_ranks > N / 2 + 1) FAIL("num_ranks=%d > N/2+1=%d\n", num_ranks, N / 2 + 1);

    int M = kernel_choose_pad_M(N, compat);
    if (kp.found[6] && kp.pad_M > 0) M = (int)kp.pad_M;
    if (pad_M_env > 0) M = pad_M_env;
    if (M < N || (M % 4) != 0 || ((M - N) % 2) != 0) {
        FAIL("padded size M=%d invalid: need M >= N=%d, M %% 4 == 0, (M - N) even\n", M, N);
    }
    if (2 * M < 3 * N && world_rank == 0) {
        fprintf(stderr, "[Kernel_MPI] WARNING: M=%d < 3N/2=%d: the quadratic product will alias\n", M, 3 * N / 2);
    }

    // f_NL selection: kernel_MPI.py overwrites kernel_output block by block, so the last
    // nonzero of (local, equil, orthoCMB, orthoLSS) wins
    FnlType fnl = FNL_NONE;
    int num_nonzero = 0;
    for (int i = 0; i < 4; i++) {
        if (kp.fnl[i] != 0.0) { fnl = (FnlType)i; num_nonzero++; }
    }
    if (fnl == FNL_NONE) FAIL("all of ZD_f_NL_local/equil/orthoCMB/orthoLSS are zero (or missing)\n");
    if (num_nonzero > 1 && world_rank == 0) {
        fprintf(stderr, "[Kernel_MPI] WARNING: %d f_NL types are nonzero; using %s (the last one, as kernel_MPI.py does)\n",
                num_nonzero, fnl_type_name(fnl));
    }

    KernelTerm terms[KERNEL_MAX_TERMS];
    int num_terms = kernel_terms_for(fnl, terms);
    if (only_term != NULL && only_term[0] != '\0') {
        int id = kernel_shape_from_name(only_term);
        if (id < 0) FAIL("KERNEL_ONLY_TERM=%s is not a kernel shape name\n", only_term);
        terms[0] = (KernelTerm){(KernelShapeId)id, 1.0};
        num_terms = 1;
    }

    // ==================================================================================
    // Topology (same Cartesian grid as Zeldovich_MPI) + layout
    // ==================================================================================
    if (zd_topology_init_standalone(num_ranks, num_z_ranks) != 0) FAIL("topology init failed\n");
    const int rank = zd_cart_rank;
    const int size_x = zd_cart_size_x, size_z = zd_cart_size_z;

    KernelLayout K;
    kernel_layout_init(&K, N, M, num_ranks, rank, size_x, size_z);

    if (rank == 0) {
        printf("========================================================================\n");
        printf("Kernel_MPI: f_NL %s (%d term%s)%s\n", fnl_type_name(fnl), num_terms, num_terms == 1 ? "" : "s",
               (only_term && only_term[0]) ? " [KERNEL_ONLY_TERM]" : "");
        printf("  N=%d  M=%d (padding %.3fx)  BoxSize=%g\n", N, M, (double)M / N, kp.box_size);
        printf("  ranks=%d (size_x=%d x size_z=%d), OMP threads/rank=%d\n",
               num_ranks, size_x, size_z, omp_get_max_threads());
        printf("  batches=%d, compat=%d, zero_nyquist=%d\n", K.num_batches, compat, zero_nyquist);
        printf("  in:  %s\n  out: %s\n", bin_dir, out_dir);
        for (int t = 0; t < num_terms; t++) {
            printf("  term %d: %s  coef=% .17g\n", t, kernel_shape(terms[t].shape)->name, terms[t].coef);
        }
        printf("========================================================================\n");
        fflush(stdout);
    }

    // ==================================================================================
    // Read this rank's primary y-slices (mirrors are built on the fly by Hermitian symmetry)
    // ==================================================================================
    const int64_t slice_elems = (int64_t)N * N;
    fftw_complex_t **phi = (fftw_complex_t**)calloc(K.my_num_pairs > 0 ? K.my_num_pairs : 1, sizeof(*phi));
    fftw_complex_t **out = (fftw_complex_t**)calloc(K.my_num_pairs > 0 ? K.my_num_pairs : 1, sizeof(*out));
    for (int b = 0; b < K.my_num_pairs; b++) {
        int y_p = kernel_batch_y_p(&K, rank, b);
        char path[4096];
        snprintf(path, sizeof(path), "%s/Phi_gaussian_rank%d_y%d.bin", bin_dir, rank, y_p);
        phi[b] = kernel_alloc_complex(slice_elems, "phi slice");
        int err = kernel_read_slice(path, N, phi[b]);
        if (err != 0) {
            fprintf(stderr, "[Kernel_MPI] ERROR: rank %d cannot read %s (%s); was the dump made with %d ranks?\n",
                    rank, path, err == 1 ? "missing" : "wrong size", num_ranks);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        if (y_p == 0 || y_p == N / 2) kernel_fix_self_conjugate_plane(N, phi[b]);
        if (zero_nyquist) kernel_zero_nyquist(N, y_p, phi[b]);

        out[b] = kernel_alloc_complex(slice_elems, "output slice");
        memset(out[b], 0, sizeof(fftw_complex_t) * (size_t)slice_elems);
    }

    // ==================================================================================
    // Kernel terms
    // ==================================================================================
    KernelPass P;
    kernel_pass_init(&P, &K, zd_comm_2d, kp.box_size, phi, out);

    if (rank == 0) {
        double gb = 1.0e-9 * sizeof(fftw_complex_t);
        printf("[Memory, rank 0] recv_buffer %.3f GB, q_half %.3f GB, phi+out %.3f GB\n",
               K.recv_total * gb, (double)(N / 2 + 1) * K.my_pencils * gb,
               2.0 * K.my_num_pairs * slice_elems * gb);
        fflush(stdout);
    }

    double t_terms = MPI_Wtime();
    for (int t = 0; t < num_terms; t++) {
        const KernelShape *shape = kernel_shape(terms[t].shape);
        double t0 = MPI_Wtime();
        run_kernel_term(&P, shape, terms[t].coef);
        MPI_Barrier(zd_comm_2d);
        if (rank == 0) {
            printf("[Kernel_MPI] term %d/%d (%s) done in %.3f s\n", t + 1, num_terms, shape->name, MPI_Wtime() - t0);
            fflush(stdout);
        }
    }
    t_terms = MPI_Wtime() - t_terms;

    // ==================================================================================
    // Write output (rank r writes exactly the y_primary slices it read)
    // ==================================================================================
    if (rank == 0) {
        if (mkdir(out_dir, 0755) != 0 && errno != EEXIST) {
            fprintf(stderr, "[Kernel_MPI] ERROR: cannot create %s\n", out_dir);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    MPI_Barrier(zd_comm_2d);
    for (int b = 0; b < K.my_num_pairs; b++) {
        int y_p = kernel_batch_y_p(&K, rank, b);
        char path[4096];
        snprintf(path, sizeof(path), "%s/phi_output_rank%d_y%d.bin", out_dir, rank, y_p);
        if (kernel_write_slice(path, N, out[b]) != 0) {
            fprintf(stderr, "[Kernel_MPI] ERROR: rank %d cannot write %s\n", rank, path);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    // ==================================================================================
    // Timing / memory summary (max over ranks)
    // ==================================================================================
    double local_t[7] = {P.t_fill_fft2d, P.t_pack, P.t_comm_fwd, P.t_zslab, P.t_comm_ret,
                         P.t_fft2d_accum, peak_rss_gb()};
    double max_t[7];
    MPI_Reduce(local_t, max_t, 7, MPI_DOUBLE, MPI_MAX, 0, zd_comm_2d);
    if (rank == 0) {
        printf("[Timing, max over ranks] fill+2D bwd %.3f s | pack %.3f s | Alltoallv fwd %.3f s | "
               "z-slab 1D+multiply %.3f s | Alltoallv ret %.3f s | 2D fwd+accumulate %.3f s\n",
               max_t[0], max_t[1], max_t[2], max_t[3], max_t[4], max_t[5]);
        printf("[Kernel_MPI] kernel terms %.3f s, total %.3f s, peak RSS (max rank) %.3f GB\n",
               t_terms, MPI_Wtime() - t_start, max_t[6]);
        fflush(stdout);
    }

    kernel_pass_free(&P);
    for (int b = 0; b < K.my_num_pairs; b++) { free(phi[b]); free(out[b]); }
    free(phi);
    free(out);
    kernel_layout_free(&K);
    zd_topology_cleanup();
    MPI_Finalize();
    return 0;
}

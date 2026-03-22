/*
 * Toy: OMP scaling test for z-loop over a single shared buffer (mimics local_y_slices).
 *
 * Replicates:
 *   - Same buffer layout (2 slices: primary + conjugate, narray*N*N each)
 *   - No omp-parallelization in Y
 *   - Parallelization over z-rows only (#pragma omp parallel for schedule(static))
 *   - Each thread: alloc RNG buffer once, get_rng_copy(ps, global_y, buf), advance to z-start
 *   - PTimerWall for "RNG setup" (alloc+copy+advance) and "Z-loop compute"
 *
 * With PLT: use param file with qPLT=1 and PLT_filename set; toy loads eig_vecs and
 * computes F,G,H from plt_get_eigenmode (same shared read as production). narray becomes 4.
 *
 * Example: OMP_NUM_THREADS=4 ./toy_zloop_omp_scaling 1024 2 param.par 10
 *
 * MPI mode (mimic production multi-node): run with mpiexec so multiple ranks each run
 * the same z-loop workload (no communication). Use to test if MPI+OMP causes poor scaling.
 *   mpiexec -n 9 ./toy_zloop_omp_scaling 1024 2 param.par 10
 *
 * Env (optional, to try to fix anti-scaling under mpiexec):
 *   TOY_ZLOOP_MAX_THREADS=N  — cap threads in the z-loop (e.g. 4) while OMP_NUM_THREADS can be higher.
 *   TOY_FIRST_TOUCH=1       — parallel first-touch of buffer so NUMA pages are faulted by the thread that will use them.
 */

#include "../src/precision.h"
#include "../src/config.h"
#include "../src/PTimer.h"
#include "../src/types.h"
#include "../src/utils/zeldovich_wrapper.h"
#include "../src/utils/plt_eigenmodes.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>

#ifndef MAX_PPD
#define MAX_PPD 65536
#endif

#ifndef SPLINE_RESOLUTION
#define SPLINE_RESOLUTION 128
#endif

static inline int64_t compute_virtual_position(int z, int x, int N, int Nhalf) {
    int z_v = (z <= Nhalf) ? z : (MAX_PPD - N + z);
    int x_v = (x <= Nhalf) ? x : (MAX_PPD - N + x);
    return 2 * ((int64_t)z_v * MAX_PPD + x_v);
}

/* Same as compute_plt_eigenmode: lookup, sign-flip, normalize, k2/(k·e) scale */
static inline int toy_plt_eigenmode(int kx, int ky, int kz, int N, double k2, eigenmode* e) {
    int ikx = (kx < 0) ? N + kx : kx;
    int iky = (ky < 0) ? N + ky : ky;
    int ikz = (kz < 0) ? N + kz : kz;
    if (ikz > N / 2) ikz = N - ikz;
    if (plt_get_eigenmode(ikx, iky, ikz, (int64_t)N, e) != 0) return 0;
    if (kz < 0) e->vec[2] = -e->vec[2];
    double e_mag = sqrt(e->vec[0]*e->vec[0] + e->vec[1]*e->vec[1] + e->vec[2]*e->vec[2]);
    if (e_mag > 0.0) {
        e->vec[0] /= e_mag; e->vec[1] /= e_mag; e->vec[2] /= e_mag;
    }
    double k_dot_e = kx*e->vec[0] + ky*e->vec[1] + kz*e->vec[2];
    double norm = (k2 > 0.0 && k_dot_e != 0.0) ? k2 / k_dot_e : 0.0;
    if (!isfinite(norm)) norm = 0.0;
    e->vec[0] *= norm; e->vec[1] *= norm; e->vec[2] *= norm;
    return 1;
}

int main(int argc, char** argv) {
    int rank = 0, size = 1;
    int use_mpi = 1;
    const char* disable_mpi_env = getenv("TOY_DISABLE_MPI");
    if (disable_mpi_env && atoi(disable_mpi_env) != 0) {
        use_mpi = 0;
    }

    if (use_mpi) {
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
    }

    if (argc < 4) {
        if (rank == 0)
            fprintf(stderr, "Usage: %s <N> <narray> <param_file> [num_y_repeats]\n", argv[0]);
        if (use_mpi) MPI_Finalize();
        return 1;
    }
    int N = atoi(argv[1]);
    int narray = atoi(argv[2]);
    const char* param_file = argv[3];
    int num_y_repeats = (argc >= 5) ? atoi(argv[4]) : 1;

    if (N <= 0 || narray <= 0 || num_y_repeats <= 0) {
        if (rank == 0)
            fprintf(stderr, "N, narray, num_y_repeats must be positive.\n");
        if (use_mpi) MPI_Finalize();
        return 1;
    }

    ParametersHandle params = zeldovich_params_create(param_file);
    if (!params) {
        if (rank == 0) fprintf(stderr, "Failed to load params from %s\n", param_file);
        if (use_mpi) MPI_Finalize();
        return 1;
    }

    double fundamental = zeldovich_params_get_fundamental(params);
    double k_cutoff = zeldovich_params_get_k_cutoff(params);
    int CornerModes = zeldovich_params_get_CornerModes(params);
    int Nhalf = N / 2;
    double Nhalf_dbl = (double)Nhalf;
    double k2_cutoff = (Nhalf_dbl * Nhalf_dbl) / (k_cutoff * k_cutoff);

    PowerSpectrumHandle ps = zeldovich_ps_create(SPLINE_RESOLUTION, params);
    if (!ps) {
        if (rank == 0) fprintf(stderr, "Failed to create PowerSpectrum\n");
        zeldovich_params_destroy(params);
        if (use_mpi) MPI_Finalize();
        return 1;
    }
    double powerlaw_index = zeldovich_params_get_Pk_powerlaw_index(params);
    if (powerlaw_index == 1000.0) {
        if (rank == 0) fprintf(stderr, "ZD_Pk_powerlaw_index not in param file\n");
        zeldovich_ps_destroy(ps);
        zeldovich_params_destroy(params);
        if (use_mpi) MPI_Finalize();
        return 1;
    }
    if (zeldovich_ps_init_powerlaw(ps, powerlaw_index, params) != 0) {
        if (rank == 0) fprintf(stderr, "Failed to init power spectrum\n");
        zeldovich_ps_destroy(ps);
        zeldovich_params_destroy(params);
        MPI_Finalize();
        return 1;
    }

    int use_plt = 0;
    double f_cluster = 0.0;
    if (zeldovich_params_get_qPLT(params)) {
        const char* PLT_filename = zeldovich_params_get_PLT_filename(params);
        if (PLT_filename && PLT_filename[0] && plt_load_eigenmodes(PLT_filename) == 0) {
            use_plt = 1;
            narray = 4;
            f_cluster = zeldovich_params_get_f_cluster(params);
        }
    }

    int64_t slice_elems = (int64_t)narray * (int64_t)N * (int64_t)N;
    size_t total_bytes = 2 * (size_t)slice_elems * sizeof(fftw_complex_t);
    fftw_complex_t* local_y_slices = NULL;
    if (posix_memalign((void**)&local_y_slices, 128, total_bytes) != 0) {
        if (rank == 0) fprintf(stderr, "Failed to allocate local_y_slices (%zu bytes)\n", total_bytes);
        zeldovich_ps_destroy(ps);
        zeldovich_params_destroy(params);
        MPI_Finalize();
        return 1;
    }

    fftw_complex_t* primary_slices   = &local_y_slices[0];
    fftw_complex_t* conjugate_slices = &local_y_slices[slice_elems];

    static PTimerWall pt_rng_setup(1);
    static PTimerWall pt_zloop(1);
    pt_rng_setup.Clear();
    pt_zloop.Clear();

    int max_t = omp_get_max_threads();
    int zloop_threads = max_t;
    const char* cap_env = getenv("TOY_ZLOOP_MAX_THREADS");
    if (cap_env) {
        int cap = atoi(cap_env);
        if (cap > 0 && cap < zloop_threads) zloop_threads = cap;
    }
    size_t rng_size = zeldovich_ps_rng_buffer_size();
    void** rng_bufs = (void**)calloc((size_t)max_t, sizeof(void*));
    int zloop_started = 0;  /* function scope so shared in parallel for */
    if (!rng_bufs) {
        if (rank == 0) fprintf(stderr, "Failed to allocate rng_bufs\n");
        free(local_y_slices);
        zeldovich_ps_destroy(ps);
        zeldovich_params_destroy(params);
        MPI_Finalize();
        return 1;
    }

    /* idea: parallel first-touch so each z-row is faulted on the thread that will write it (NUMA) */
    if (getenv("TOY_FIRST_TOUCH") && atoi(getenv("TOY_FIRST_TOUCH")) != 0) {
#pragma omp parallel for schedule(static) num_threads(zloop_threads)
        for (int z = 0; z < N; z++) {
            for (int a = 0; a < narray; a++) {
                int64_t p0 = (int64_t)N * ((int64_t)z + (int64_t)N * (int64_t)a);
                primary_slices[p0][0] = 0;
                primary_slices[p0][1] = 0;
            }
            int z_mirror = (z == 0) ? 0 : N - z;
            for (int a = 0; a < narray; a++) {
                int64_t c0 = (int64_t)N * ((int64_t)z_mirror + (int64_t)N * (int64_t)a);
                conjugate_slices[c0][0] = 0;
                conjugate_slices[c0][1] = 0;
            }
        }
    }

    if (rank == 0) {
        printf("Toy z-loop OMP scaling: N=%d narray=%d num_y_repeats=%d OMP_threads=%d zloop_threads=%d use_plt=%d MPI_ranks=%d\n",
               N, narray, num_y_repeats, max_t, zloop_threads, use_plt, size);
        printf("  Buffer: %zu bytes (same layout as production local_y_slices)\n", total_bytes);
        if (getenv("TOY_FIRST_TOUCH")) printf("  TOY_FIRST_TOUCH=%s\n", getenv("TOY_FIRST_TOUCH"));
        fflush(stdout);
    }

    for (int y_rep = 0; y_rep < num_y_repeats; y_rep++) {
        int global_y = 1;
        memset(local_y_slices, 0, total_bytes);

        pt_rng_setup.Start(0);
        zloop_started = 0;

        void* local_rng_buf = NULL;
        int64_t nskip = 0;
        int rng_ready = 0;
#pragma omp parallel for schedule(static) num_threads(zloop_threads) firstprivate(local_rng_buf, nskip, rng_ready)
        for (int z = 0; z < N; z++) {
                if (!rng_ready) {
                    int tid = omp_get_thread_num();
                    local_rng_buf = malloc(rng_size);
                    rng_bufs[tid] = local_rng_buf;
                    zeldovich_ps_get_rng_copy(ps, global_y, local_rng_buf);
                    int64_t vstart = compute_virtual_position(z, 0, N, Nhalf) / 2;
                    if (vstart > 0)
                        zeldovich_ps_advance_rng_buffer(local_rng_buf, vstart);
                    rng_ready = 1;

                    int prev;
#pragma omp atomic capture
                    prev = zloop_started++;
                    if (prev == 0) {
                        pt_rng_setup.Stop(0);
                        pt_zloop.Start(0);
                    }
                }

                if (z == Nhalf + 1 && N < MAX_PPD) {
                    int64_t skip_amount = (int64_t)(MAX_PPD - N) * (int64_t)MAX_PPD;
                    nskip += skip_amount;
                }

                int ky = (global_y > Nhalf) ? global_y - N : global_y;
                int kz = (z > Nhalf) ? z - N : z;
                int abs_ky = (ky < 0) ? -ky : ky;
                int abs_kz = (kz < 0) ? -kz : kz;

                for (int x = 0; x < N; x++) {
                    if (x == Nhalf + 1 && N < MAX_PPD) {
                        nskip += (int64_t)(MAX_PPD - N);
                    }

                    int kx = (x > Nhalf) ? x - N : x;
                    int k2_int = kx * kx + ky * ky + kz * kz;
                    double k2 = (double)k2_int;
                    int abs_kx = (kx < 0) ? -kx : kx;
                    int is_nyquist = (abs_kx == Nhalf || abs_ky == Nhalf || abs_kz == Nhalf);

                    double D_re = 0.0, D_im = 0.0;
                    if ((k2 == 0.0) || is_nyquist || (!CornerModes && (double)k2_int >= k2_cutoff)) {
                        nskip++;
                    } else {
                        double k2_phys = k2 * fundamental * fundamental;
                        double kmag = sqrt(k2_phys);
                        if (nskip > 0) {
                            zeldovich_ps_advance_rng_buffer(local_rng_buf, nskip);
                            nskip = 0;
                        }
                        zeldovich_ps_cgauss_from_buffer(local_rng_buf, ps, kmag, &D_re, &D_im);
                    }

                    double F_re = 0.0, F_im = 0.0, G_re = 0.0, G_im = 0.0, H_re = 0.0, H_im = 0.0;
                    double f_vel = 1.0;
                    if (D_re != 0.0 || D_im != 0.0) {
                        if (k2 == 0.0) f_vel = 0.0;
                        else {
                            double factor = 1.0 / (k2 * fundamental);
                            if (use_plt) {
                                eigenmode e;
                                if (toy_plt_eigenmode(kx, ky, kz, N, k2, &e)) {
                                    f_vel = (sqrt(1. + 24. * e.val * f_cluster) - 1.) * 0.25;
                                    F_re = -e.vec[0] * factor * D_im; F_im =  e.vec[0] * factor * D_re;
                                    G_re = -e.vec[1] * factor * D_im; G_im =  e.vec[1] * factor * D_re;
                                    H_re = -e.vec[2] * factor * D_im; H_im =  e.vec[2] * factor * D_re;
                                } else {
                                    F_re = -kx * factor * D_im; F_im = kx * factor * D_re;
                                    G_re = -ky * factor * D_im; G_im = ky * factor * D_re;
                                    H_re = -kz * factor * D_im; H_im = kz * factor * D_re;
                                }
                            } else {
                                F_re = -kx * factor * D_im; F_im = kx * factor * D_re;
                                G_re = -ky * factor * D_im; G_im = ky * factor * D_re;
                                H_re = -kz * factor * D_im; H_im = kz * factor * D_re;
                            }
                        }
                    }

                    int x_mirror = (x == 0) ? 0 : N - x;
                    int z_mirror = (z == 0) ? 0 : N - z;
                    /* Store like production store_prim_conj: arrays 0,1 always; 2,3 if narray>=4 */
                    int64_t p0 = (int64_t)x + (int64_t)N * ((int64_t)z + (int64_t)N * 0);
                    int64_t p1 = (int64_t)x + (int64_t)N * ((int64_t)z + (int64_t)N * 1);
                    int64_t c0 = (int64_t)x_mirror + (int64_t)N * ((int64_t)z_mirror + (int64_t)N * 0);
                    int64_t c1 = (int64_t)x_mirror + (int64_t)N * ((int64_t)z_mirror + (int64_t)N * 1);
                    primary_slices[p0][0] = (real_t)(D_re - F_im); primary_slices[p0][1] = (real_t)(D_im + F_re);
                    primary_slices[p1][0] = (real_t)(G_re - H_im); primary_slices[p1][1] = (real_t)(G_im + H_re);
                    conjugate_slices[c0][0] = (real_t)(D_re + F_im); conjugate_slices[c0][1] = (real_t)(F_re - D_im);
                    conjugate_slices[c1][0] = (real_t)(G_re + H_im); conjugate_slices[c1][1] = (real_t)(H_re - G_im);
                    if (narray >= 4) {
                        int64_t p2 = (int64_t)x + (int64_t)N * ((int64_t)z + (int64_t)N * 2);
                        int64_t p3 = (int64_t)x + (int64_t)N * ((int64_t)z + (int64_t)N * 3);
                        int64_t c2 = (int64_t)x_mirror + (int64_t)N * ((int64_t)z_mirror + (int64_t)N * 2);
                        int64_t c3 = (int64_t)x_mirror + (int64_t)N * ((int64_t)z_mirror + (int64_t)N * 3);
                        primary_slices[p2][0] = (real_t)(-F_im * f_vel); primary_slices[p2][1] = (real_t)(F_re * f_vel);
                        primary_slices[p3][0] = (real_t)((G_re - H_im) * f_vel); primary_slices[p3][1] = (real_t)((G_im + H_re) * f_vel);
                        conjugate_slices[c2][0] = (real_t)(F_im * f_vel); conjugate_slices[c2][1] = (real_t)(F_re * f_vel);
                        conjugate_slices[c3][0] = (real_t)((G_re + H_im) * f_vel); conjugate_slices[c3][1] = (real_t)((H_re - G_im) * f_vel);
                    }
                }
            }

        pt_zloop.Stop(0);

        for (int t = 0; t < max_t; t++) {
            if (rng_bufs[t]) {
                free(rng_bufs[t]);
                rng_bufs[t] = NULL;
            }
        }
    }

    free(rng_bufs);

    double rng_s = pt_rng_setup.Elapsed();
    double zloop_s = pt_zloop.Elapsed();
    double total = rng_s + zloop_s;
    printf("[Rank %d] [PTimerWall] %d Y-slice repeats: RNG_setup=%.6fs Z-loop=%.6fs Sum=%.6fs\n",
           rank, num_y_repeats, rng_s, zloop_s, total);
    fflush(stdout);

    free(local_y_slices);
    if (use_plt) plt_free_eigenmodes();
    zeldovich_ps_destroy(ps);
    zeldovich_params_destroy(params);
    if (use_mpi) MPI_Finalize();
    return 0;
}

/*
 * Minimal MPI toy to debug FFTW wisdom + the same batched 2D / 1D plan path as production
 * (src/fft/fft_setup.c), with wisdom handling modeled on zeldovich-PLT Setup_FFTW:
 *   - Fixed filename in cwd: fftw_toy.wisdom (like fftw_zeldovich.wisdom there)
 *   - fftw_import_wisdom_from_filename, print return to stderr
 *   - create plans (setup_fftw_plans_full without USE_FFTW_WISDOM — wisdom done here)
 *   - fftw_export_wisdom_to_filename, print return to stderr
 *
 * MPI: rank 0 reads the file; wisdom string is broadcast so all ranks share identical
 * planner state before threading + plan creation (zeldovich is single-process only).
 *
 * Build:  make toy-wisdom-debug
 *
 * Run from the directory where fftw_toy.wisdom should live (e.g. toy/bin_files).
 */

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <mpi.h>
#include <omp.h>

extern "C" {
#include "config.h"
#include "fft/fft_setup.h"
#include "precision.h"
}

namespace fs = std::filesystem;

/* Same pattern as zeldovich.cpp Setup_FFTW, plus MPI broadcast of wisdom string. */
static void toy_wisdom_import_mpi(int rank, MPI_Comm comm)
{
    const fs::path wisdom_file("fftw_toy.wisdom");
    char *wisdom_str = nullptr;
    size_t wisdom_len = 0;

    if (rank == 0) {
        int wisdom_exists = FFTW_IMPORT_WISDOM_FROM_FILENAME(wisdom_file.c_str());
        fprintf(stderr,
                "FFTW Wisdom import from file \"%s\" returned %d (%s).\n",
                wisdom_file.string().c_str(),
                wisdom_exists,
                wisdom_exists == 1 ? "success" : "failure");
        fflush(stderr);
        if (!wisdom_exists) {
            FFTW_FORGET_WISDOM();
        }
        wisdom_str = FFTW_EXPORT_WISDOM_TO_STRING();
        if (wisdom_str != nullptr) {
            wisdom_len = std::strlen(wisdom_str);
        }
    }

    MPI_Bcast(&wisdom_len, 1, MPI_UNSIGNED_LONG_LONG, 0, comm);

    if (wisdom_len > 0) {
        if (rank != 0) {
            wisdom_str = (char *)std::malloc(wisdom_len + 1);
            if (wisdom_str == nullptr) {
                fprintf(stderr, "[toy-wisdom] rank %d: malloc failed for wisdom broadcast\n", rank);
                MPI_Abort(comm, 1);
            }
        }
        MPI_Bcast(wisdom_str, (int)wisdom_len, MPI_CHAR, 0, comm);
        wisdom_str[wisdom_len] = '\0';
        FFTW_FORGET_WISDOM();
        FFTW_IMPORT_WISDOM_FROM_STRING(wisdom_str);
        if (rank == 0) {
            FFTW_FREE(wisdom_str);
        } else {
            std::free(wisdom_str);
        }
    }
}

static void toy_wisdom_export_rank0(int rank)
{
    if (rank != 0) {
        return;
    }
    const fs::path wisdom_file("fftw_toy.wisdom");
    int ret = FFTW_EXPORT_WISDOM_TO_FILENAME(wisdom_file.c_str());
    fprintf(stderr, "FFTW Wisdom export to file \"%s\" returned %d.\n", wisdom_file.string().c_str(), ret);
    fflush(stderr);
}

int main(int argc, char** argv)
{
    int provided = 0;
    const int required = MPI_THREAD_FUNNELED;
    if (MPI_Init_thread(nullptr, nullptr, required, &provided) != MPI_SUCCESS) {
        fprintf(stderr, "MPI_Init_thread failed\n");
        return 1;
    }
    if (provided < required) {
        fprintf(stderr, "Need MPI_THREAD_FUNNELED for FFTW+OpenMP\n");
        MPI_Finalize();
        return 1;
    }

    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc < 4) {
        if (rank == 0)
            fprintf(stderr,
                    "Usage: %s <N> <narray> <num_y_repeats>\n"
                    "  Run from cwd where fftw_toy.wisdom should appear (e.g. bin_files).\n"
                    "  Example: mpiexec -n 64 %s 2048 4 17\n",
                    argv[0], argv[0]);
        MPI_Finalize();
        return 1;
    }

    int N = std::atoi(argv[1]);
    int narray = std::atoi(argv[2]);
    int num_y_repeats = std::atoi(argv[3]);

    if (N <= 0 || (N % 4) != 0 || narray <= 0 || num_y_repeats <= 0) {
        if (rank == 0)
            fprintf(stderr, "Invalid args: need N>0, N%%4==0, narray>0, num_y_repeats>0\n");
        MPI_Finalize();
        return 1;
    }

    const int64_t slice_elems = (int64_t)narray * (int64_t)N * (int64_t)N;
    const size_t total_elems = (size_t)(2 * slice_elems);
    const size_t total_bytes = total_elems * sizeof(fftw_complex_t);

    fftw_complex_t* local_y_slices = nullptr;
    if (posix_memalign((void**)&local_y_slices, (size_t)ALIGN_BYTES, total_bytes) != 0) {
        if (rank == 0)
            fprintf(stderr, "posix_memalign failed (%zu bytes)\n", total_bytes);
        MPI_Finalize();
        return 1;
    }

    fftw_complex_t* primary_slices = &local_y_slices[0];
    fftw_complex_t* conjugate_slices = &local_y_slices[slice_elems];

    for (size_t i = 0; i < total_elems; i++) {
        local_y_slices[i][0] = (float)(i % 97);
        local_y_slices[i][1] = (float)((i * 7) % 91);
    }

    fftw_plan_t plan_2d{};
    fftw_plan_t plan_1d_y{};

    if (rank == 0) {
        printf("[toy-wisdom] N=%d narray=%d repeats=%d ranks=%d OMP=%d\n",
               N, narray, num_y_repeats, size, omp_get_max_threads());
        fflush(stdout);
    }

    toy_wisdom_import_mpi(rank, MPI_COMM_WORLD);
    setup_fftw_plans_full(N, narray, primary_slices, &plan_2d, &plan_1d_y);
    toy_wisdom_export_rank0(rank);

    double t0 = MPI_Wtime();
    for (int yrep = 0; yrep < num_y_repeats; yrep++) {
        (void)yrep;
        FFTW_EXECUTE_DFT(plan_2d, primary_slices, primary_slices);
        FFTW_EXECUTE_DFT(plan_2d, conjugate_slices, conjugate_slices);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    double t1 = MPI_Wtime();

    if (rank == 0) {
        printf("[toy-wisdom] 2D execute loop: %.6f s (per rank, %d repeats x 2 slabs)\n",
               t1 - t0, num_y_repeats);
        fflush(stdout);
    }

    FFTW_DESTROY_PLAN(plan_2d);
    FFTW_DESTROY_PLAN(plan_1d_y);

#ifdef USE_DOUBLE_PRECISION
    fftw_cleanup_threads();
#else
    fftwf_cleanup_threads();
#endif

    std::free(local_y_slices);

    MPI_Finalize();
    return 0;
}

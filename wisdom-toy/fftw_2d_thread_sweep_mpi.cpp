/*
every rank runs the same batched 2D FFT
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <mpi.h>
#include <fftw3.h>
#include <omp.h>
#include <unistd.h>

#ifndef ALIGN_BYTES
#define ALIGN_BYTES 4096
#endif

int main(int argc, char **argv)
{
    int provided = 0;
    MPI_Init_thread(NULL, NULL, MPI_THREAD_FUNNELED, &provided);

    int rank = 0, nranks = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nranks);

    char host[256];
    host[0] = '\0';
    if (gethostname(host, sizeof(host) - 1) != 0) {
        std::snprintf(host, sizeof(host), "unknown");
    }

    const int N = std::atoi(argv[1]);
    const int narray = std::atoi(argv[2]);
    const int fftw_threads = std::atoi(argv[3]);
    const int repeats = (argc >= 5) ? std::atoi(argv[4]) : 20;

    const unsigned flags = FFTW_MEASURE;

    omp_set_num_threads(fftw_threads);

    if (fftwf_init_threads() == 0) {
        if (rank == 0) {
            std::fprintf(stderr, "fftwf_init_threads failed\n");
        }
        MPI_Finalize();
        return 1;
    }

    const char *wisdom_file = std::getenv("FFTW_WISDOM_FILE");
    if (wisdom_file != nullptr && wisdom_file[0] != '\0') {
        const int wok = fftwf_import_wisdom_from_filename(wisdom_file);
        if (rank == 0) {
            std::fprintf(stderr,
                         "[FFTW-WISDOM] import \"%s\": %d (%s)\n",
                         wisdom_file,
                         wok,
                         wok ? "ok" : "missing or parse failed (will plan from scratch)");
            std::fflush(stderr);
        }
        if (!wok) {
            fftwf_forget_wisdom();
        }
    }

    fftwf_plan_with_nthreads(fftw_threads);

    const size_t plane_elems = static_cast<size_t>(N) * static_cast<size_t>(N);
    const size_t total_elems = plane_elems * static_cast<size_t>(narray);
    const size_t nbytes = total_elems * sizeof(fftwf_complex);

    fftwf_complex *buf = static_cast<fftwf_complex*>(fftwf_malloc(nbytes));
    if (buf == nullptr) {
        if (rank == 0) {
            std::fprintf(stderr, "fftwf_malloc failed\n");
        }
        fftwf_cleanup_threads();
        MPI_Finalize();
        return 1;
    }
    std::memset(buf, 0, nbytes);
    for (int a = 0; a < narray; ++a) {
        buf[static_cast<size_t>(a) * plane_elems][0] = 1.0f;
    }

    int n[2] = { N, N };
    const int idist = static_cast<int>(plane_elems);

    MPI_Barrier(MPI_COMM_WORLD);
    const double t_plan0 = omp_get_wtime();
    fftwf_plan plan = fftwf_plan_many_dft(
        2,
        n,
        narray,
        buf,
        nullptr,
        1,
        idist,
        buf,
        nullptr,
        1,
        idist,
        FFTW_FORWARD,
        flags);
    const double t_plan1 = omp_get_wtime();
    MPI_Barrier(MPI_COMM_WORLD);

    if (!plan) {
        if (rank == 0) {
            std::fprintf(stderr, "fftwf_plan_many_dft failed\n");
        }
        fftwf_free(buf);
        fftwf_cleanup_threads();
        MPI_Finalize();
        return 1;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    const double t_exec0 = omp_get_wtime();
    for (int r = 0; r < repeats; ++r) {
        fftwf_execute_dft(plan, buf, buf);
    }
    const double t_exec1 = omp_get_wtime();

    const double exec_s = t_exec1 - t_exec0;
    const double per_call_ms = 1000.0 * exec_s / static_cast<double>(repeats);

    const double gb = static_cast<double>(nbytes) / 1.0e9;
    const double eff_gbs = (2.0 * gb * static_cast<double>(repeats)) / exec_s;

    std::printf(
        "rank=%d/%d host=%s N=%d narray=%d fftw_threads=%d repeats=%d flags=%s\n",
        rank,
        nranks,
        host,
        N,
        narray,
        fftw_threads,
        repeats,
        (flags & FFTW_MEASURE) ? "MEASURE" : "ESTIMATE");
    std::printf("rank=%d plan_time_s=%.6f execute_total_s=%.6f execute_per_call_ms=%.6f effective_GB_s=%.3f\n",
                rank,
                t_plan1 - t_plan0,
                exec_s,
                per_call_ms,
                eff_gbs);
    std::fflush(stdout);

    if (wisdom_file != nullptr && wisdom_file[0] != '\0' && rank == 0) {
        const int exok = fftwf_export_wisdom_to_filename(wisdom_file);
        std::fprintf(stderr,
                     "[FFTW-WISDOM] export \"%s\": %d (%s)\n",
                     wisdom_file,
                     exok,
                     exok ? "ok" : "failed");
        std::fflush(stderr);
    }

    fftwf_destroy_plan(plan);
    fftwf_free(buf);
    fftwf_cleanup_threads();
    MPI_Finalize();
    return 0;
}

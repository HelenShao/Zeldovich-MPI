/*
 * batched 2D complex FFT (same layout as fft_setup.c /
 * generate_hermitian_slice_pair_local), threaded FFTW, sweep OpenMP thread counts.
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fftw3.h>
#include <omp.h>

#ifndef ALIGN_BYTES
#define ALIGN_BYTES 4096
#endif

int main(int argc, char **argv)
{
    const int N = std::atoi(argv[1]);
    const int narray = std::atoi(argv[2]);
    const int fftw_threads = std::atoi(argv[3]);
    const int repeats = (argc >= 5) ? std::atoi(argv[4]) : 20;
    const unsigned flags = FFTW_PATIENT;

    omp_set_num_threads(fftw_threads);
    if (fftwf_init_threads() == 0) {
        std::fprintf(stderr, "fftwf_init_threads failed\n");
        return 1;
    }

    const char *wisdom_file = std::getenv("FFTW_WISDOM_FILE");
    if (wisdom_file != nullptr && wisdom_file[0] != '\0') {
        const int wok = fftwf_import_wisdom_from_filename(wisdom_file);
        std::fprintf(stderr,
                     "[FFTW-WISDOM] import \"%s\": %d (%s)\n",
                     wisdom_file,
                     wok,
                     wok ? "ok" : "missing or parse failed (will plan from scratch)");
        std::fflush(stderr);
    }

    fftwf_plan_with_nthreads(fftw_threads);

    const size_t plane_elems = static_cast<size_t>(N) * static_cast<size_t>(N);
    const size_t total_elems = plane_elems * static_cast<size_t>(narray);
    const size_t nbytes = total_elems * sizeof(fftwf_complex);


    //fftwf_complex *buf = nullptr;
    // if (posix_memalign(reinterpret_cast<void **>(&buf), ALIGN_BYTES, nbytes) != 0) {
    //     std::fprintf(stderr, "posix_memalign failed\n");
    //     fftwf_cleanup_threads();
    //     return 1;
    // }
    fftwf_complex *buf = static_cast<fftwf_complex*>(fftwf_malloc(nbytes));
    if (buf == nullptr) {
        std::fprintf(stderr, "fftwf_malloc failed\n");
        fftwf_cleanup_threads();
        return 1;
    }
    std::memset(buf, 0, nbytes);
    /* first touch - core that touches first will control the memory location, to be NUMA aware, want the thread that does the fft to do the first touch */
    for (int a = 0; a < narray; ++a) {
        buf[static_cast<size_t>(a) * plane_elems][0] = 1.0f;
    }

    int n[2] = { N, N };
    const int idist = static_cast<int>(plane_elems);

    const double t_plan0 = omp_get_wtime();
    fftwf_plan plan = fftwf_plan_many_dft(
        2,            // 2D transform
        n,            // n[2]: transform size = {N, N}
        narray,       // # of 2D transforms in batch
        buf,          // input base pointer
        nullptr,      // inembed: use n[] as physical dims (contiguous planes)
        1,            // istride: adjacent elements in fastest dim are contiguous
        idist,        // idist: element offset between successive planes
        buf,          // output base pointer (in-place)
        nullptr,      // onembed: use n[] as physical dims for output
        1,            // ostride: contiguous output in fastest dim
        idist,        // odist: element offset between output planes
        FFTW_FORWARD, 
        flags);       // planner flags: FFTW_PATIENT
    const double t_plan1 = omp_get_wtime();

    if (!plan) {
        std::fprintf(stderr, "fftwf_plan_many_dft failed\n");
        fftwf_free(buf);
        fftwf_cleanup_threads();
        return 1;
    }

    const double t_exec0 = omp_get_wtime();
    for (int r = 0; r < repeats; ++r) {
        fftwf_execute_dft(plan, buf, buf);
    }
    const double t_exec1 = omp_get_wtime();

    const double exec_s = t_exec1 - t_exec0;
    const double per_call_ms = 1000.0 * exec_s / static_cast<double>(repeats);

    /* Rough "effective GB/s" for complex array: read+write in-place ~ 2 * nbytes per exec */
    const double gb = static_cast<double>(nbytes) / 1.0e9;
    const double eff_gbs = (2.0 * gb * static_cast<double>(repeats)) / exec_s;

    const char *flag_name = "UNKNOWN";
    if (flags & FFTW_PATIENT) {
        flag_name = "PATIENT";
    } else if (flags & FFTW_MEASURE) {
        flag_name = "MEASURE";
    } else if (flags & FFTW_ESTIMATE) {
        flag_name = "ESTIMATE";
    }

    std::printf(
        "N=%d narray=%d fftw_threads=%d omp_max=%d repeats=%d flags=%s\n",
        N,
        narray,
        fftw_threads,
        omp_get_max_threads(),
        repeats,
        flag_name);
    std::printf("plan_time_s=%.6f\n", t_plan1 - t_plan0);
    std::printf("execute_total_s=%.6f execute_per_call_ms=%.6f effective_GB_s=%.3f\n",
                exec_s,
                per_call_ms,
                eff_gbs);

    if (wisdom_file != nullptr && wisdom_file[0] != '\0') {
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
    return 0;
}

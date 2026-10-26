/*
 * Staged v3 batched 2D complex FFT benchmark:
 * - one reusable full-plane stage buffer
 * - copy in/out uses 64-byte row chunks (8 complex values)
 * - full 2D FFT runs on fully packed stage buffer
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <fftw3.h>
#include <omp.h>

#ifndef ALIGN_BYTES
#define ALIGN_BYTES 4096
#endif

namespace {
constexpr int kChunkComplex = 8; // 8 complex values * 8B = 64B

void copy_plane_chunked(fftwf_complex *dst, const fftwf_complex *src, int N)
{
    #pragma omp parallel for schedule(static)
    for (int y = 0; y < N; ++y) {
        fftwf_complex *dst_row = dst + static_cast<size_t>(y) * static_cast<size_t>(N);
        const fftwf_complex *src_row = src + static_cast<size_t>(y) * static_cast<size_t>(N);
        for (int x = 0; x < N; x += kChunkComplex) {
            // step in 8 complex numbers at a time (cache line size)
            // paste into dst_row + x, x = 0, 8, 16, ...
            // then go to next row (y)
            const int len = (x + kChunkComplex <= N) ? kChunkComplex : (N - x);
            std::memcpy(dst_row + x, src_row + x, static_cast<size_t>(len) * sizeof(fftwf_complex));
        }
    }
}
} 

int main(int argc, char **argv)
{
    if (argc < 4) {
        std::fprintf(stderr,
                     "Usage: %s N NARRAY FFTW_THREADS [REPEATS] [PATIENT|MEASURE|ESTIMATE]\n",
                     argv[0]);
        return 2;
    }

    const int N = std::atoi(argv[1]);
    const int narray = std::atoi(argv[2]);
    const int fftw_threads = std::atoi(argv[3]);
    const int repeats = (argc >= 5) ? std::atoi(argv[4]) : 20;
    const char *planner_arg = (argc >= 6) ? argv[5] : "PATIENT";

    if (N <= 0 || narray <= 0 || fftw_threads <= 0 || repeats <= 0) {
        std::fprintf(stderr, "N, NARRAY, FFTW_THREADS, and REPEATS must all be > 0\n");
        return 2;
    }

    unsigned flags = FFTW_PATIENT;
    if (std::strcmp(planner_arg, "PATIENT") == 0) {
        flags = FFTW_PATIENT;
    } else if (std::strcmp(planner_arg, "MEASURE") == 0) {
        flags = FFTW_MEASURE;
    } else if (std::strcmp(planner_arg, "ESTIMATE") == 0) {
        flags = FFTW_ESTIMATE;
    } else {
        std::fprintf(stderr,
                     "Invalid planner flag: %s (expected PATIENT, MEASURE, or ESTIMATE)\n",
                     planner_arg);
        return 2;
    }

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
    const size_t plane_bytes = plane_elems * sizeof(fftwf_complex);
    const size_t nbytes = total_elems * sizeof(fftwf_complex);

    fftwf_complex *buf = static_cast<fftwf_complex *>(fftwf_malloc(nbytes));
    if (buf == nullptr) {
        std::fprintf(stderr, "fftwf_malloc failed\n");
        fftwf_cleanup_threads();
        return 1;
    }

    // Parallel first touch so pages are mapped by worker threads
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < total_elems; ++i) {
        buf[i][0] = 0.0f;
        buf[i][1] = 0.0f;
    }

    fftwf_complex *stage_buf = nullptr;
    {
        void *raw = nullptr;
        if (posix_memalign(&raw, ALIGN_BYTES, plane_bytes) != 0) {
            std::fprintf(stderr, "posix_memalign failed for stage buffer\n");
            fftwf_free(buf);
            fftwf_cleanup_threads();
            return 1;
        }
        stage_buf = static_cast<fftwf_complex *>(raw);
        std::memset(stage_buf, 0, plane_bytes);
    }

    // fftw plan creation
    const double t_plan0 = omp_get_wtime();
    fftwf_plan plan = fftwf_plan_dft_2d(
        N,
        N,
        stage_buf,
        stage_buf,
        FFTW_FORWARD,
        flags);
    const double t_plan1 = omp_get_wtime();

    // fftw execution 
    double copy_in_s = 0.0;
    double copy_out_s = 0.0;
    double fft_only_s = 0.0;

    const double t_exec0 = omp_get_wtime();
    for (int r = 0; r < repeats; ++r) {
        for (int a = 0; a < narray; ++a) {
            fftwf_complex *plane = &buf[static_cast<size_t>(a) * plane_elems];
            const double t_copy_in0 = omp_get_wtime();
            copy_plane_chunked(stage_buf, plane, N);
            const double t_copy_in1 = omp_get_wtime();
            copy_in_s += (t_copy_in1 - t_copy_in0);

            const double t_fft0 = omp_get_wtime();
            fftwf_execute_dft(plan, stage_buf, stage_buf);
            const double t_fft1 = omp_get_wtime();
            fft_only_s += (t_fft1 - t_fft0);

            const double t_copy_out0 = omp_get_wtime();
            copy_plane_chunked(plane, stage_buf, N);
            const double t_copy_out1 = omp_get_wtime();
            copy_out_s += (t_copy_out1 - t_copy_out0);
        }
    }
    const double t_exec1 = omp_get_wtime();

    const double exec_s = t_exec1 - t_exec0;
    const double per_call_ms = 1000.0 * exec_s / static_cast<double>(repeats);
    const double slices = static_cast<double>(repeats) * static_cast<double>(narray);
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
    std::printf("mode=STAGED_V3_STRIDED64 align_bytes=%d chunk_complex=%d\n",
                ALIGN_BYTES,
                kChunkComplex);
    std::printf("plan_time_s=%.6f\n", t_plan1 - t_plan0);
    std::printf("execute_total_s=%.6f execute_per_call_ms=%.6f effective_GB_s=%.3f\n",
                exec_s,
                per_call_ms,
                eff_gbs);
    std::printf("copy_in_total_s=%.6f copy_out_total_s=%.6f fft_execute_only_total_s=%.6f\n",
                copy_in_s,
                copy_out_s,
                fft_only_s);
    std::printf(
        "copy_in_per_repeat_ms=%.6f copy_out_per_repeat_ms=%.6f fft_execute_only_per_repeat_ms=%.6f "
        "copy_in_per_slice_ms=%.6f copy_out_per_slice_ms=%.6f fft_execute_only_per_slice_ms=%.6f\n",
        1000.0 * copy_in_s / static_cast<double>(repeats),
        1000.0 * copy_out_s / static_cast<double>(repeats),
        1000.0 * fft_only_s / static_cast<double>(repeats),
        1000.0 * copy_in_s / slices,
        1000.0 * copy_out_s / slices,
        1000.0 * fft_only_s / slices);

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
    std::free(stage_buf);
    fftwf_free(buf);
    fftwf_cleanup_threads();
    return 0;
}

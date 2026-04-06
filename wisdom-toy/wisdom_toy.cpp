/*
 * Minimal FFTW wisdom import/export test (single precision, same as parent Makefile default).
 *
 * 1) Batched 2D DFTs via FFTW_PLAN_MANY_DFT: narray identical N×N transforms (matches fft_setup.c).
 * 2) 1D complex DFT via plan_and_execute_1d()
 *
 * Wisdom: fftwf uses one in-memory wisdom database per process. A single import before any
 * planning and a single export after all plans both apply to the 2D batch and the 1D transform;
 * the same file can hold entries for all float plans you create (FFTW merges on export).
 *
 * Default wisdom file: ./wisdom_toy.wisdom
 */

#include <cstdio>
#include <cstdlib>

#include <fftw3.h>

#ifndef FFTW_PLAN_MANY_DFT
#define FFTW_PLAN_MANY_DFT(rank, n, howmany, in, inembed, istride, idist, out, onembed, ostride, odist, sign, flags) \
    fftwf_plan_many_dft(rank, n, howmany, in, inembed, istride, idist, out, onembed, ostride, odist, sign, flags)
#endif
#ifndef FFTW_PLAN_DFT_1D
#define FFTW_PLAN_DFT_1D(n, in, out, sign, flags) fftwf_plan_dft_1d(n, in, out, sign, flags)
#endif
#ifndef PLAN_FFT_1D
#define PLAN_FFT_1D(n, in, out, sign, flags) FFTW_PLAN_DFT_1D(n, in, out, sign, flags)
#endif
#ifndef FFTW_IMPORT_WISDOM_FROM_FILENAME
#define FFTW_IMPORT_WISDOM_FROM_FILENAME(fn) fftwf_import_wisdom_from_filename(fn)
#endif
#ifndef FFTW_EXPORT_WISDOM_TO_FILENAME
#define FFTW_EXPORT_WISDOM_TO_FILENAME(fn) fftwf_export_wisdom_to_filename(fn)
#endif
#ifndef FFTW_FORGET_WISDOM
#define FFTW_FORGET_WISDOM() fftwf_forget_wisdom()
#endif
#ifndef FFTW_ALLOC_COMPLEX
#define FFTW_ALLOC_COMPLEX(n) fftwf_alloc_complex(n)
#endif
#ifndef FFTW_FREE
#define FFTW_FREE(p) fftwf_free(p)
#endif
#ifndef FFTW_EXECUTE_DFT
#define FFTW_EXECUTE_DFT(p, in, out) fftwf_execute_dft(p, in, out)
#endif
#ifndef FFTW_DESTROY_PLAN
#define FFTW_DESTROY_PLAN(p) fftwf_destroy_plan(p)
#endif

#ifndef WISDOM_FILE_DEFAULT
#define WISDOM_FILE_DEFAULT "wisdom_toy.wisdom"
#endif

/* Start of plane `a` in a contiguous [narray][N][N] layout (C row-major within each plane). */
static inline fftwf_complex *plane_base(fftwf_complex *buf, int a, int N)
{
    return buf + static_cast<size_t>(a) * static_cast<size_t>(N) * static_cast<size_t>(N);
}

/*
 * Plan (FFTW_MEASURE) and execute one out-of-place length-n complex 1D DFT.
 * Uses the same fftwf in-memory wisdom as any prior fftwf_plan_* in this process.
 * Returns 1 on success, 0 if planning failed.
 */
static int plan_and_execute_1d(int n, fftwf_complex *in, fftwf_complex *out)
{
    fftwf_plan p = PLAN_FFT_1D(n, in, out, FFTW_FORWARD, FFTW_MEASURE);
    if (!p) {
        std::fprintf(stderr, "PLAN_FFT_1D failed (n=%d)\n", n);
        return 0;
    }
    FFTW_EXECUTE_DFT(p, in, out);
    FFTW_DESTROY_PLAN(p);
    return 1;
}

int main(int argc, char **argv)
{
    const char *wisdom_path = (argc >= 2) ? argv[1] : WISDOM_FILE_DEFAULT;

    const int wisdom_import = FFTW_IMPORT_WISDOM_FROM_FILENAME(wisdom_path);
    std::fprintf(stderr,
                 "FFTW wisdom import from \"%s\": %d (%s)\n",
                 wisdom_path,
                 wisdom_import,
                 wisdom_import ? "ok" : "failed or missing");
    std::fflush(stderr);

    if (!wisdom_import) {
        FFTW_FORGET_WISDOM();
    }

    /* Same roles as fft_setup.c: narray batched N×N 2D transforms. */
    const int N = 512;
    const int narray = 4;
    const int n1d = 65536;

    const size_t elems_per_plane = static_cast<size_t>(N) * static_cast<size_t>(N);
    const size_t buf_2d_elems = elems_per_plane * static_cast<size_t>(narray);

    fftwf_complex *buf_2d = FFTW_ALLOC_COMPLEX(buf_2d_elems);
    fftwf_complex *in1 = FFTW_ALLOC_COMPLEX(static_cast<size_t>(n1d));
    fftwf_complex *out1 = FFTW_ALLOC_COMPLEX(static_cast<size_t>(n1d));

    if (!buf_2d || !in1 || !out1) {
        std::fprintf(stderr, "allocation failed\n");
        return 1;
    }

    /* Initialize each 2D plane separately in the stacked buffer. */
    for (int a = 0; a < narray; ++a) {
        fftwf_complex *p = plane_base(buf_2d, a, N);
        for (size_t i = 0; i < elems_per_plane; ++i) {
            p[i][0] = 0.0f;
            p[i][1] = 0.0f;
        }
        /* DC bin scaled by (a+1) for easy identification. */
        p[0][0] = static_cast<float>(a + 1);
    }

    for (int i = 0; i < n1d; ++i) {
        in1[i][0] = static_cast<float>(i % 17);
        in1[i][1] = 0.0f;
    }

    /*
     * FFTW_PLAN_MANY_DFT (same pattern as src/fft/fft_setup.c):
     *   rank=2, n[] = { N, N }  — each transform is N×N; n[1] is the contiguous (stride-1) dimension.
     *   howmany = narray       — number of such 2D transforms.
     *   istride=1, idist=N*N   — within a plane, consecutive along last index; next plane starts +N*N.
     *   inembed/onembed NULL   — use n[] as full storage dimensions for each slab.
     * In-place: same buffer for in and out (matches production plan_buffer usage).
     */
    int n[2] = { N, N };
    fftwf_plan plan_many_2d = FFTW_PLAN_MANY_DFT(
        2,
        n,
        narray,
        buf_2d, nullptr, 1, static_cast<int>(elems_per_plane),
        buf_2d, nullptr, 1, static_cast<int>(elems_per_plane),
        FFTW_FORWARD,
        FFTW_MEASURE);

    if (!plan_many_2d) {
        std::fprintf(stderr, "FFTW_PLAN_MANY_DFT (2D batched) failed\n");
        FFTW_FREE(buf_2d);
        FFTW_FREE(in1);
        FFTW_FREE(out1);
        return 1;
    }
    FFTW_EXECUTE_DFT(plan_many_2d, buf_2d, buf_2d);
    FFTW_DESTROY_PLAN(plan_many_2d);

    if (!plan_and_execute_1d(n1d, in1, out1)) {
        FFTW_FREE(buf_2d);
        FFTW_FREE(in1);
        FFTW_FREE(out1);
        return 1;
    }

    const int wisdom_export = FFTW_EXPORT_WISDOM_TO_FILENAME(wisdom_path);
    std::fprintf(stderr,
                 "FFTW wisdom export to \"%s\": %d (%s)\n",
                 wisdom_path,
                 wisdom_export,
                 wisdom_export ? "ok" : "failed");
    std::fflush(stderr);

    FFTW_FREE(buf_2d);
    FFTW_FREE(in1);
    FFTW_FREE(out1);

    return wisdom_export ? 0 : 2;
}

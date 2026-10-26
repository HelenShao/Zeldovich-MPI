// ====================================================================================
// KERNEL FILE I/O + Y-SLICE HERMITIAN HELPERS
// ====================================================================================

#include "kernel_io.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

int kernel_detect_N(const char *path)
{
    struct stat st;
    if (stat(path, &st) != 0) return -1;
    long long n_complex = (long long)st.st_size / (2 * (long long)sizeof(double));
    if (n_complex * 2 * (long long)sizeof(double) != (long long)st.st_size) return -1;
    long long N = llround(sqrt((double)n_complex));
    if (N * N != n_complex) return -1;
    return (int)N;
}

#ifndef USE_DOUBLE_PRECISION
// Single precision: files stay float64, so convert through a bounded double buffer
#define KERNEL_IO_CHUNK ((size_t)1 << 16)   // complex values per chunk (1 MiB of doubles)
#endif

int kernel_read_slice(const char *path, int N, fftw_complex_t *slice)
{
    FILE *fp = fopen(path, "rb");
    if (fp == NULL) return 1;
    size_t nvals = (size_t)N * (size_t)N;
#ifdef USE_DOUBLE_PRECISION
    // fftw_complex (double precision) is layout-compatible with interleaved doubles
    size_t nread = fread(slice, 2 * sizeof(double), nvals, fp);
#else
    double *buf = (double*)malloc(KERNEL_IO_CHUNK * 2 * sizeof(double));
    size_t nread = 0;
    while (buf != NULL && nread < nvals) {
        size_t want = (nvals - nread < KERNEL_IO_CHUNK) ? nvals - nread : KERNEL_IO_CHUNK;
        size_t got = fread(buf, 2 * sizeof(double), want, fp);
        for (size_t i = 0; i < got; i++) {
            slice[nread + i][0] = (real_t)buf[2 * i];
            slice[nread + i][1] = (real_t)buf[2 * i + 1];
        }
        nread += got;
        if (got != want) break;
    }
    free(buf);
#endif
    int extra = fgetc(fp);
    fclose(fp);
    return (nread == nvals && extra == EOF) ? 0 : 2;
}

int kernel_write_slice(const char *path, int N, const fftw_complex_t *slice)
{
    FILE *fp = fopen(path, "wb");
    if (fp == NULL) return 1;
    size_t nvals = (size_t)N * (size_t)N;
#ifdef USE_DOUBLE_PRECISION
    size_t nwritten = fwrite(slice, 2 * sizeof(double), nvals, fp);
#else
    double *buf = (double*)malloc(KERNEL_IO_CHUNK * 2 * sizeof(double));
    size_t nwritten = 0;
    while (buf != NULL && nwritten < nvals) {
        size_t want = (nvals - nwritten < KERNEL_IO_CHUNK) ? nvals - nwritten : KERNEL_IO_CHUNK;
        for (size_t i = 0; i < want; i++) {
            buf[2 * i]     = (double)slice[nwritten + i][0];
            buf[2 * i + 1] = (double)slice[nwritten + i][1];
        }
        size_t put = fwrite(buf, 2 * sizeof(double), want, fp);
        nwritten += put;
        if (put != want) break;
    }
    free(buf);
#endif
    int close_err = fclose(fp);
    return (nwritten == nvals && close_err == 0) ? 0 : 2;
}

// Element (x, z) of a [z][x] slice
#define SL(x, z) slice[(size_t)(z) * (size_t)N + (size_t)(x)]

void kernel_fix_self_conjugate_plane(int N, fftw_complex_t *slice)
{
    int Nh = N / 2;

    // Mirror the lower half onto the upper half: P[(N-x)%N, N-z] = conj(P[x, z])
    // for x in [0, N), z in [1, N/2). Reads rows z < N/2, writes rows z > N/2.
    // Rows z = 0 and z = N/2 map onto themselves and are handled below.
    for (int z = 1; z < Nh; z++) {
        for (int x = 0; x < N; x++) {
            int xm = (N - x) % N;
            SL(xm, N - z)[0] =  SL(x, z)[0];
            SL(xm, N - z)[1] = -SL(x, z)[1];
        }
    }
    // P[N - x, 0] = conj(P[x, 0]) and P[N - x, N/2] = conj(P[x, N/2]) for x in [1, N/2)
    for (int x = 1; x < Nh; x++) {
        SL(N - x, 0)[0]  =  SL(x, 0)[0];
        SL(N - x, 0)[1]  = -SL(x, 0)[1];
        SL(N - x, Nh)[0] =  SL(x, Nh)[0];
        SL(N - x, Nh)[1] = -SL(x, Nh)[1];
    }
    // Force real at DC/Nyquist corners
    SL(0, 0)[1] = 0.0;
    SL(0, Nh)[1] = 0.0;
    SL(Nh, 0)[1] = 0.0;
    SL(Nh, Nh)[1] = 0.0;
}

void kernel_zero_nyquist(int N, int y_p, fftw_complex_t *slice)
{
    int Nh = N / 2;
    if (y_p == Nh) {
        memset(slice, 0, sizeof(fftw_complex_t) * (size_t)N * (size_t)N);
        return;
    }
    for (int i = 0; i < N; i++) {
        SL(Nh, i)[0] = SL(Nh, i)[1] = 0.0;   // kx = N/2
        SL(i, Nh)[0] = SL(i, Nh)[1] = 0.0;   // kz = N/2
    }
}

#undef SL

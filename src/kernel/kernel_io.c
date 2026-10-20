// ====================================================================================
// KERNEL FILE I/O + Y-SLICE HERMITIAN HELPERS
// ====================================================================================

#include "kernel_io.h"
#include <dirent.h>
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

int kernel_max_dumped_rank(const char *dir)
{
    DIR *d = opendir(dir);
    if (d == NULL) return -1;
    int max_rank = -1;
    struct dirent *ent;
    while ((ent = readdir(d)) != NULL) {
        int r = -1, y = -1;
        char tail[8] = {0};
        if (sscanf(ent->d_name, "Phi_gaussian_rank%d_y%d%7s", &r, &y, tail) == 3 &&
            strcmp(tail, ".bin") == 0 && r > max_rank) {
            max_rank = r;
        }
    }
    closedir(d);
    return max_rank;
}

int kernel_read_slice(const char *path, int N, fftw_complex_t *slice)
{
    FILE *fp = fopen(path, "rb");
    if (fp == NULL) return 1;
    size_t nvals = (size_t)N * (size_t)N;
    // fftw_complex (double precision) is layout-compatible with interleaved doubles
    size_t nread = fread(slice, 2 * sizeof(double), nvals, fp);
    int extra = fgetc(fp);
    fclose(fp);
    return (nread == nvals && extra == EOF) ? 0 : 2;
}

int kernel_write_slice(const char *path, int N, const fftw_complex_t *slice)
{
    FILE *fp = fopen(path, "wb");
    if (fp == NULL) return 1;
    size_t nvals = (size_t)N * (size_t)N;
    size_t nwritten = fwrite(slice, 2 * sizeof(double), nvals, fp);
    int close_err = fclose(fp);
    return (nwritten == nvals && close_err == 0) ? 0 : 2;
}

// Element (x, z) of a [z][x] slice
#define SL(x, z) slice[(size_t)(z) * (size_t)N + (size_t)(x)]

void kernel_fix_self_conjugate_plane(int N, fftw_complex_t *slice)
{
    int Nh = N / 2;

    // P[(-x)%N, N - z] = conj(P[x, z]) for all x, z in [1, N/2)   (reads z < N/2, writes z > N/2)
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

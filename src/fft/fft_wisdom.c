#include "fft_wisdom.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void build_rank_wisdom_filename(int rank, char *out, size_t out_size)
{
    snprintf(out, out_size, "%s_rank%04d", FFTW_WISDOM_FILENAME, rank);
}

void fft_wisdom_import_per_rank(int rank, MPI_Comm comm)
{
    (void)comm; // outdated
    int imported = 0;
    char filename[256];
    char *wisdom_str = NULL;
    size_t wisdom_len = 0;

    build_rank_wisdom_filename(rank, filename, sizeof(filename));
    imported = FFTW_IMPORT_WISDOM_FROM_FILENAME(filename);

    if (!imported) {
        printf("[FFTW-WISDOM] Rank %d: no usable wisdom file '%s'; planning from scratch.\n",
               rank, filename);
        fflush(stdout);
        return;
    }

    wisdom_str = FFTW_EXPORT_WISDOM_TO_STRING();
    wisdom_len = wisdom_str ? strlen(wisdom_str) : 0;
    if (wisdom_str != NULL) {
        FFTW_FREE(wisdom_str);
    }

    // if (wisdom_len <= 128) {
    //     printf("[FFTW-WISDOM] Rank %d: '%s' has no plan entries (len=%zu); "
    //            "discarding and planning from scratch.\n",
    //            rank, filename, wisdom_len);
    //     fflush(stdout);
    //     FFTW_FORGET_WISDOM();
    //     return;
    // }

    printf("[FFTW-WISDOM] Rank %d: imported wisdom from '%s' (%s precision, len=%zu).\n",
           rank, filename, PRECISION_NAME, wisdom_len);
    fflush(stdout);

    // MPI_Bcast(&imported, 1, MPI_INT, 0, comm);
    // MPI_Bcast(&wisdom_len, 1, MPI_UNSIGNED_LONG_LONG, 0, comm);
    // MPI_Bcast(wisdom_str, (int)wisdom_len, MPI_CHAR, 0, comm);
}

void fft_wisdom_export_per_rank(int rank)
{
    char filename[256];
    char *ws = FFTW_EXPORT_WISDOM_TO_STRING();
    if (ws) {
        printf("[FFTW-WISDOM] Rank %d: accumulated wisdom = %zu bytes\n",
               rank, strlen(ws));
        printf("[FFTW-WISDOM] Rank %d: wisdom string = %s\n",
               rank, ws);
        fflush(stdout);
        FFTW_FREE(ws);
    }

    build_rank_wisdom_filename(rank, filename, sizeof(filename));
    int ok = FFTW_EXPORT_WISDOM_TO_FILENAME(filename);
    if (!ok) {
        fprintf(stderr,
                "[FFTW-WISDOM] Rank %d: failed to export wisdom to '%s' (%s precision)\n",
                rank, filename, PRECISION_NAME);
        fflush(stderr);
    } else {
        printf("[FFTW-WISDOM] Rank %d: exported wisdom to '%s' (%s precision)\n",
               rank, filename, PRECISION_NAME);
        fflush(stdout);
    }
}


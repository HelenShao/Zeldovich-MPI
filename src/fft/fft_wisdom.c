#include "fft_wisdom.h"

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int fft_wisdom_import_rank0_broadcast_local(int rank, MPI_Comm comm, const char *local_wisdom_dir)
{
    (void)local_wisdom_dir;

    char *wisdom_str = NULL;
    size_t wisdom_len = 0;
    int wisdom_from_fftw_alloc = 0;

    if (rank == 0) {
        const int imported = FFTW_IMPORT_WISDOM_FROM_FILENAME(FFTW_WISDOM_FILENAME);
        if (!imported) {
            fprintf(stderr, "[FFTW-WISDOM] Rank 0: failed to import wisdom from '%s' (%s precision)\n",
                    FFTW_WISDOM_FILENAME, PRECISION_NAME);
            fflush(stderr);
            return -1;
        }

        wisdom_str = FFTW_EXPORT_WISDOM_TO_STRING();
        if (wisdom_str == NULL) {
            fprintf(stderr, "[FFTW-WISDOM] Rank 0: FFTW_EXPORT_WISDOM_TO_STRING returned NULL\n");
            fflush(stderr);
            return -1;
        }
        wisdom_from_fftw_alloc = 1;
        wisdom_len = strlen(wisdom_str);
        if (wisdom_len == 0u) {
            fprintf(stderr, "[FFTW-WISDOM] Rank 0: exported wisdom string is empty\n");
            fflush(stderr);
            FFTW_FREE(wisdom_str);
            return -1;
        }
    }

    MPI_Bcast(&wisdom_len, 1, MPI_UNSIGNED_LONG_LONG, 0, comm);
    if (wisdom_len == 0u) {
        if (rank == 0 && wisdom_str != NULL && wisdom_from_fftw_alloc) {
            FFTW_FREE(wisdom_str);
        }
        fprintf(stderr, "[FFTW-WISDOM] Rank %d: received empty wisdom\n", rank);
        fflush(stderr);
        return -1;
    }
    if (wisdom_len > (size_t)INT_MAX) {
        if (rank == 0 && wisdom_str != NULL && wisdom_from_fftw_alloc) {
            FFTW_FREE(wisdom_str);
        }
        fprintf(stderr, "[FFTW-WISDOM] Rank %d: wisdom too large for MPI_Bcast (%zu)\n", rank,
                wisdom_len);
        fflush(stderr);
        return -1;
    }

    if (rank != 0) {
        wisdom_str = (char *)malloc(wisdom_len + 1u);
        if (wisdom_str == NULL) {
            fprintf(stderr, "[FFTW-WISDOM] Rank %d: failed to allocate wisdom buffer (%zu bytes)\n", rank,
                    wisdom_len + 1u);
            fflush(stderr);
            return -1;
        }
    }

    MPI_Bcast(wisdom_str, (int)wisdom_len, MPI_CHAR, 0, comm);
    wisdom_str[wisdom_len] = '\0';

    FFTW_FORGET_WISDOM();
    const int imported = FFTW_IMPORT_WISDOM_FROM_STRING(wisdom_str);
    if (!imported) {
        fprintf(stderr, "[FFTW-WISDOM] Rank %d: failed to import broadcast wisdom string (%s precision)\n",
                rank, PRECISION_NAME);
        fflush(stderr);
        if (rank == 0 && wisdom_from_fftw_alloc) {
            FFTW_FREE(wisdom_str);
        } else {
            free(wisdom_str);
        }
        return -1;
    }

    if (rank == 0 && wisdom_from_fftw_alloc) {
        FFTW_FREE(wisdom_str);
    } else {
        free(wisdom_str);
    }
    return 0;
}

#include "fft_wisdom.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef FFTW_WISDOM_PREVIEW_MAX
#define FFTW_WISDOM_PREVIEW_MAX 512u
#endif

#ifndef FFTW_WISDOM_PATH_BUFSZ
#define FFTW_WISDOM_PATH_BUFSZ 4096
#endif

static const char *fftw_wisdom_file_path(char *buf)
{
    const char *env = getenv("FFTW_WISDOM_FILE");
    if (env != NULL && env[0] != '\0') {
        if (snprintf(buf, FFTW_WISDOM_PATH_BUFSZ, "%s", env) >= FFTW_WISDOM_PATH_BUFSZ) {
            fprintf(stderr, "[FFTW-WISDOM] FFTW_WISDOM_FILE path too long; using '%s'.\n",
                    FFTW_WISDOM_FILENAME);
            fflush(stderr);
            return FFTW_WISDOM_FILENAME;
        }
        return buf;
    }
    return FFTW_WISDOM_FILENAME;
}

void fft_wisdom_import_broadcast(int rank, MPI_Comm comm)
{
    (void)comm; // outdated
    int imported = 0;
    char filename[256];
    char *wisdom_str = NULL;
    size_t wisdom_len = 0;

    if (rank == 0) {
        char pathbuf[FFTW_WISDOM_PATH_BUFSZ];
        const char *wpath = fftw_wisdom_file_path(pathbuf);

        imported = FFTW_IMPORT_WISDOM_FROM_FILENAME(wpath);

        if (!imported) {
            // Ensure we start from a clean state on all ranks.
            FFTW_FORGET_WISDOM();
            printf("[FFTW-WISDOM] No existing '%s' for %s precision; "
                   "planning from scratch, will create it later.\n",
                   wpath, PRECISION_NAME);
            fflush(stdout);
        } else {
            printf("[FFTW-WISDOM] Imported wisdom from '%s' (%s precision).\n",
                   wpath, PRECISION_NAME);
            fflush(stdout);
        }

        // Export current (possibly empty) wisdom state to a string to share
        // exactly the same wisdom with all ranks.
        wisdom_str = FFTW_EXPORT_WISDOM_TO_STRING();
        if (wisdom_str != NULL) {
            wisdom_len = strlen(wisdom_str);
        } else {
            wisdom_len = 0;
        }
    }

    // Broadcast the length first.
    MPI_Bcast(&wisdom_len, 1, MPI_UNSIGNED_LONG_LONG, 0, comm);

    if (wisdom_len > 0) {
        if (rank != 0) {
            wisdom_str = (char *)malloc(wisdom_len + 1);
            if (!wisdom_str) {
                fprintf(stderr,
                        "[FFTW-WISDOM] Rank %d: failed to allocate buffer for wisdom broadcast\n",
                        rank);
                MPI_Abort(comm, 1);
            }
        }

        MPI_Bcast(wisdom_str, (int)wisdom_len, MPI_CHAR, 0, comm);
        wisdom_str[wisdom_len] = '\0';

        // All ranks (including 0) import the shared wisdom string.
        FFTW_FORGET_WISDOM();
        FFTW_IMPORT_WISDOM_FROM_STRING(wisdom_str);

        if (rank == 0) {
            printf("[FFTW-WISDOM] Broadcasted wisdom to all ranks (len=%zu).\n",
                   wisdom_len);
            fflush(stdout);
        }
    } else {
        // No wisdom available; all ranks already have a clean state.
        if (rank == 0) {
            printf("[FFTW-WISDOM] No wisdom to broadcast; all ranks start fresh.\n");
            fflush(stdout);
        }
    }

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

    char *pre = FFTW_EXPORT_WISDOM_TO_STRING();
    if (pre != NULL) {
        size_t n = strlen(pre);
        printf("[FFTW-WISDOM] Rank 0: pre-export wisdom length %zu bytes (%s precision)\n",
               n, PRECISION_NAME);
        if (n <= FFTW_WISDOM_PREVIEW_MAX) {
            printf("[FFTW-WISDOM] Rank 0: wisdom (full):\n%.*s\n", (int)n, pre);
        } else {
            printf("[FFTW-WISDOM] Rank 0: wisdom (first %u bytes):\n%.*s\n... (%zu more bytes)\n",
                   (unsigned)FFTW_WISDOM_PREVIEW_MAX, (int)FFTW_WISDOM_PREVIEW_MAX, pre,
                   n - (size_t)FFTW_WISDOM_PREVIEW_MAX);
        }
        fflush(stdout);
        FFTW_FREE(pre);
    } else {
        fprintf(stderr,
                "[FFTW-WISDOM] Rank 0: FFTW_EXPORT_WISDOM_TO_STRING returned NULL (%s precision)\n",
                PRECISION_NAME);
        fflush(stderr);
    }

    char pathbuf[FFTW_WISDOM_PATH_BUFSZ];
    const char *wpath = fftw_wisdom_file_path(pathbuf);

    int ok = FFTW_EXPORT_WISDOM_TO_FILENAME(wpath);
    if (!ok) {
        fprintf(stderr,
                "[FFTW-WISDOM] Rank 0: failed to export wisdom to '%s' (%s precision)\n",
                wpath, PRECISION_NAME);
        fflush(stderr);
    } else {
        printf("[FFTW-WISDOM] Rank 0: exported wisdom to '%s' (%s precision)\n",
               wpath, PRECISION_NAME);
        fflush(stdout);
    }
}


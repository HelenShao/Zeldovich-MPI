#include "fft_wisdom.h"

#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>

#define DEFAULT_LOCAL_WISDOM_DIR "/dev/shm/Abacus_wisdom"

#ifndef FFTW_WISDOM_PREVIEW_MAX
#define FFTW_WISDOM_PREVIEW_MAX 512u // cap on printed wisdom string (during debug)
#endif

int fft_wisdom_import_rank0_broadcast_local(int rank, MPI_Comm comm, const char *local_wisdom_dir)
{
    const char *target_dir = local_wisdom_dir;
    if (target_dir == NULL || target_dir[0] == '\0') {
        target_dir = DEFAULT_LOCAL_WISDOM_DIR;
    }

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

    // broadcast length of wisdom str for mem allocation on each rank 
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

    // broadcast wisdom string for each rank
    MPI_Bcast(wisdom_str, (int)wisdom_len, MPI_CHAR, 0, comm);
    wisdom_str[wisdom_len] = '\0';

    // create wisdom dir on /dev/shm if it doesn't exist (0775 = permissins)
    if (mkdir(target_dir, 0775) != 0 && errno != EEXIST) {
        fprintf(stderr, "[FFTW-WISDOM] Rank %d: failed to create local wisdom dir '%s': %s\n", rank,
                target_dir, strerror(errno));
        fflush(stderr);
        if (rank == 0 && wisdom_from_fftw_alloc) {
            FFTW_FREE(wisdom_str);
        } else {
            free(wisdom_str);
        }
        return -1;
    }

    char local_path[PATH_MAX];
    const int n = snprintf(local_path, sizeof(local_path), "%s/fftw_wisdom_rank_%d.wisdom", target_dir, rank);
    if (n <= 0 || (size_t)n >= sizeof(local_path)) {
        fprintf(stderr, "[FFTW-WISDOM] Rank %d: local wisdom path too long for dir '%s'\n", rank, target_dir);
        fflush(stderr);
        if (rank == 0 && wisdom_from_fftw_alloc) {
            FFTW_FREE(wisdom_str);
        } else {
            free(wisdom_str);
        }
        return -1;
    }

    FILE *wf = fopen(local_path, "wb");
    if (wf == NULL) {
        fprintf(stderr, "[FFTW-WISDOM] Rank %d: failed to open '%s' for write: %s\n", rank, local_path,
                strerror(errno));
        fflush(stderr);
        if (rank == 0 && wisdom_from_fftw_alloc) {
            FFTW_FREE(wisdom_str);
        } else {
            free(wisdom_str);
        }
        return -1;
    }
    if (fwrite(wisdom_str, 1, wisdom_len, wf) != wisdom_len) {
        fprintf(stderr, "[FFTW-WISDOM] Rank %d: failed to write wisdom to '%s'\n", rank, local_path);
        fflush(stderr);
        fclose(wf);
        if (rank == 0 && wisdom_from_fftw_alloc) {
            FFTW_FREE(wisdom_str);
        } else {
            free(wisdom_str);
        }
        return -1;
    }
    fclose(wf);

    const int imported_local = FFTW_IMPORT_WISDOM_FROM_FILENAME(local_path);
    if (!imported_local) {
        fprintf(stderr, "[FFTW-WISDOM] Rank %d: failed to import local wisdom '%s' (%s precision)\n", rank,
                local_path, PRECISION_NAME);
        fflush(stderr);
        if (rank == 0 && wisdom_from_fftw_alloc) {
            FFTW_FREE(wisdom_str);
        } else {
            free(wisdom_str);
        }
        return -1;
    }

    // free temporary buffer on rank 0 after received bcast
    // two different frees cuz rank0 wisdom_str came from FFTW_EXPORT_WISDOM_TO_STRING()
    if (rank == 0 && wisdom_from_fftw_alloc) {
        FFTW_FREE(wisdom_str);
    } else {
        free(wisdom_str);
    }
    return 0;
}

void fft_wisdom_export_rank0(int rank)
{
    if (rank != 0) {
        return;
    }

    char *pre = FFTW_EXPORT_WISDOM_TO_STRING();
    if (pre != NULL) {
        const size_t n = strlen(pre);
        printf("[FFTW-WISDOM] Rank 0: pre-export wisdom length %zu bytes (%s precision)\n", n,
               PRECISION_NAME);
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

    const int ok = FFTW_EXPORT_WISDOM_TO_FILENAME(FFTW_WISDOM_FILENAME);
    if (!ok) {
        fprintf(stderr,
                "[FFTW-WISDOM] Rank 0: failed to export wisdom to '%s' (%s precision)\n",
                FFTW_WISDOM_FILENAME, PRECISION_NAME);
        fflush(stderr);
    } else {
        printf("[FFTW-WISDOM] Rank 0: exported wisdom to '%s' (%s precision)\n",
               FFTW_WISDOM_FILENAME, PRECISION_NAME);
        fflush(stdout);
    }
}

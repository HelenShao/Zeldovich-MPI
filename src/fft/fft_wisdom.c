#include "fft_wisdom.h"

#include <errno.h>
#include <limits.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

static char zd_wisdom_rank0_path[PATH_MAX] = FFTW_WISDOM_BASENAME;
static bool zd_wisdom_preflight_broadcast_done_flag = false;

void zd_wisdom_set_dir(const char *dir)
{
    if (dir == NULL || dir[0] == '\0') {
        snprintf(zd_wisdom_rank0_path, sizeof(zd_wisdom_rank0_path), "%s", FFTW_WISDOM_BASENAME);
        return;
    }
    const int n = snprintf(
        zd_wisdom_rank0_path,
        sizeof(zd_wisdom_rank0_path),
        "%s/%s",
        dir,
        FFTW_WISDOM_BASENAME
    );
    if (n <= 0 || (size_t)n >= sizeof(zd_wisdom_rank0_path)) {
        snprintf(zd_wisdom_rank0_path, sizeof(zd_wisdom_rank0_path), "%s", FFTW_WISDOM_BASENAME);
    }
}

const char *zd_wisdom_rank0_file(void)
{
    return zd_wisdom_rank0_path;
}

int zd_wisdom_ensure_dir_rank0(void)
{
    const char *path = zd_wisdom_rank0_file();
    const char *slash = strrchr(path, '/');
    if (slash == NULL) {
        return 0;
    }

    const size_t dir_len = (size_t)(slash - path);
    if (dir_len == 0u || dir_len >= (size_t)PATH_MAX) {
        fprintf(stderr, "[FFTW-WISDOM] Rank 0: wisdom directory path too long\n");
        fflush(stderr);
        return -1;
    }

    char dir[PATH_MAX];
    memcpy(dir, path, dir_len);
    dir[dir_len] = '\0';

    if (mkdir(dir, 0775) != 0 && errno != EEXIST) {
        fprintf(stderr, "[FFTW-WISDOM] Rank 0: failed to create wisdom dir '%s': %s\n", dir,
                strerror(errno));
        fflush(stderr);
        return -1;
    }
    return 0;
}

void zd_wisdom_set_preflight_broadcast_done(bool done)
{
    zd_wisdom_preflight_broadcast_done_flag = done;
}

bool zd_wisdom_preflight_broadcast_done(void)
{
    return zd_wisdom_preflight_broadcast_done_flag;
}

static int zd_fftw_threads_once(int rank)
{
    static int done = 0;
    if (done) {
        return 0;
    }
    if (FFTW_INIT_THREADS() == 0) {
        fprintf(stderr, "[FFTW-WISDOM] Rank %d: FFTW_INIT_THREADS failed (%s precision)\n", rank,
                PRECISION_NAME);
        fflush(stderr);
        return -1;
    }
    done = 1;
    return 0;
}

static int fft_wisdom_import_string_on_all_ranks(
    int rank,
    char *wisdom_str,
    size_t wisdom_len,
    int wisdom_from_fftw_alloc
)
{
    (void)wisdom_len;
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

/** Rank 0: use FFTW in-memory wisdom when present; else import from preflight file. */
static int zd_wisdom_rank0_export(const char *local_wisdom_dir, char **wisdom_str_out,
                                                  size_t *wisdom_len_out, int *wisdom_from_fftw_alloc_out)
{
    char *wisdom_str = FFTW_EXPORT_WISDOM_TO_STRING();
    if (wisdom_str != NULL && wisdom_str[0] != '\0') {
        *wisdom_str_out = wisdom_str;
        *wisdom_len_out = strlen(wisdom_str);
        *wisdom_from_fftw_alloc_out = 1;
        return 0;
    }
    if (wisdom_str != NULL) {
        FFTW_FREE(wisdom_str);
    }

    if (local_wisdom_dir != NULL && local_wisdom_dir[0] != '\0') {
        zd_wisdom_set_dir(local_wisdom_dir);
    }

    const char *wisdom_file = zd_wisdom_rank0_file();
    const int imported = FFTW_IMPORT_WISDOM_FROM_FILENAME(wisdom_file);
    if (!imported) {
        fprintf(stderr,
                "[FFTW-WISDOM] Rank 0: no in-memory wisdom and failed to import '%s' (%s precision)\n",
                wisdom_file, PRECISION_NAME);
        fflush(stderr);
        return -1;
    }

    wisdom_str = FFTW_EXPORT_WISDOM_TO_STRING();
    if (wisdom_str == NULL) {
        fprintf(stderr, "[FFTW-WISDOM] Rank 0: FFTW_EXPORT_WISDOM_TO_STRING returned NULL\n");
        fflush(stderr);
        return -1;
    }
    *wisdom_from_fftw_alloc_out = 1;
    *wisdom_len_out = strlen(wisdom_str);
    if (*wisdom_len_out == 0u) {
        fprintf(stderr, "[FFTW-WISDOM] Rank 0: exported wisdom string is empty\n");
        fflush(stderr);
        FFTW_FREE(wisdom_str);
        return -1;
    }
    *wisdom_str_out = wisdom_str;
    return 0;
}

int fft_wisdom_broadcast_from_rank0(int rank, MPI_Comm comm, const char *local_wisdom_dir)
{
    if (zd_wisdom_preflight_broadcast_done()) {
        return 0;
    }

    if (zd_fftw_threads_once(rank) != 0) {
        return -1;
    }

    char *wisdom_str = NULL;
    size_t wisdom_len = 0;
    int wisdom_from_fftw_alloc = 0;

    if (rank == 0) {
        if (zd_wisdom_rank0_export(local_wisdom_dir, &wisdom_str, &wisdom_len,
                                    &wisdom_from_fftw_alloc) != 0) {
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

    const int rc = fft_wisdom_import_string_on_all_ranks(rank, wisdom_str, wisdom_len, wisdom_from_fftw_alloc);
    if (rc == 0) {
        zd_wisdom_set_preflight_broadcast_done(true);
    }
    return rc;
}

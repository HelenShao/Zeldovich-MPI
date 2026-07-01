#include "ic_stage_api.h"

#include "config.h"
#include <limits.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "fft/wisdom_rank0.h"
#include "fft/fft_wisdom.h"
#include "ic_embed_flags.h"
#include "zeldovich_mpi_driver.h"
#include "zeldovich_wrapper.h"

bool zeldovich_ic_embedded = false;
ZeldovichEmbedParamHeader zeldovich_embed_param_header = {NULL, 0};
const char *zd_ic_wisdom_save_dir = NULL;

extern "C" {

void IC_InitStage(int from_abacus_host)
{
    zeldovich_ic_embedded = (from_abacus_host != 0);
    zd_ic_wisdom_save_dir = NULL;
    zd_wisdom_set_preflight_broadcast_done(false);
    zd_wisdom_set_dir(NULL);
}

static int wisdom_narray_from_params(ParametersHandle params, int *narray_out)
{
    if (!params || !narray_out) {
        return 1;
    }

    const int qdensity = zeldovich_params_get_qdensity(params);
    if (qdensity == 2) {
        *narray_out = 1;
    } else {
        const int qPLT = zeldovich_params_get_qPLT(params);
        *narray_out = qPLT ? 4 : 2;
    }
    return 0;
}

static int wisdom_preflight_from_param_buffer(
    const char *bytes,
    size_t len,
    const char *param_path,
    const char *wisdom_save_dir
)
{
    int world_rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    const int save_to_disk = (wisdom_save_dir != NULL && wisdom_save_dir[0] != '\0');
    zd_wisdom_set_preflight_broadcast_done(false);

    if (world_rank == 0) {
        if (bytes == NULL || len < 2 || param_path == NULL || param_path[0] == '\0') {
            fprintf(stderr, "wisdom_preflight_from_param_buffer: invalid header or param_path\n");
            MPI_Barrier(MPI_COMM_WORLD);
            return 1;
        }

        ParametersHandle params = zd_params_from_buffer(bytes, len, param_path);
        if (!params) {
            fprintf(stderr,
                    "wisdom_preflight_from_param_buffer: failed to parse parameters from %s\n",
                    param_path);
            MPI_Barrier(MPI_COMM_WORLD);
            return 1;
        }

        const int64_t ppd64 = zeldovich_params_get_ppd(params);
        int narray = 0;
        const int narray_rc = wisdom_narray_from_params(params, &narray);

        if (save_to_disk) {
            zd_wisdom_set_dir(wisdom_save_dir);
        } else {
            const char *parsed_wisdom_dir = zeldovich_params_get_local_wisdom_dir(params);
            if (parsed_wisdom_dir != NULL && parsed_wisdom_dir[0] != '\0') {
                zd_wisdom_set_dir(parsed_wisdom_dir);
            } else {
                zd_wisdom_set_dir(NULL);
            }
        }

        zeldovich_params_destroy(params);
        if (narray_rc != 0) {
            MPI_Barrier(MPI_COMM_WORLD);
            return 1;
        }

        if (save_to_disk && zd_wisdom_ensure_dir_rank0() != 0) {
            MPI_Barrier(MPI_COMM_WORLD);
            return 1;
        }

        if (ppd64 <= 0 || ppd64 > (int64_t)INT_MAX) {
            fprintf(stderr,
                    "wisdom_preflight_from_param_buffer: invalid ppd=%lld from %s\n",
                    (long long)ppd64,
                    param_path);
            MPI_Barrier(MPI_COMM_WORLD);
            return 1;
        }

        const int N = (int)ppd64;
        const size_t nbytes = (size_t)N * (size_t)N * sizeof(fftw_complex_t);
        fftw_complex_t *plan_buffer = NULL;
        if (posix_memalign((void **)&plan_buffer, ALIGN_BYTES, nbytes) != 0) {
            fprintf(stderr,
                    "wisdom_preflight_from_param_buffer: posix_memalign failed (%zu bytes)\n",
                    nbytes);
            MPI_Barrier(MPI_COMM_WORLD);
            return 1;
        }

        fftw_plan_t plan_2d = NULL;
        fftw_plan_t plan_1d = NULL;
        const int rc =
            wisdom_rank0_plans(N, narray, plan_buffer, &plan_2d, &plan_1d, save_to_disk);
        if (plan_2d) {
            FFTW_DESTROY_PLAN(plan_2d);
        }
        if (plan_1d) {
            FFTW_DESTROY_PLAN(plan_1d);
        }
        free(plan_buffer);

        if (rc != 0) {
            MPI_Barrier(MPI_COMM_WORLD);
            return 1;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    return 0;
}

int IC_Rank0Wisdom(const char *param_file)
{
    int world_rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (world_rank != 0) {
        MPI_Barrier(MPI_COMM_WORLD);
        return 0;
    }

    if (param_file == NULL || param_file[0] == '\0') {
        fprintf(stderr, "IC_Rank0Wisdom: empty param_file\n");
        MPI_Barrier(MPI_COMM_WORLD);
        return 1;
    }

    ParametersHandle params = zeldovich_params_create(param_file);
    if (!params) {
        fprintf(stderr, "IC_Rank0Wisdom: failed to load parameters from %s\n", param_file);
        MPI_Barrier(MPI_COMM_WORLD);
        return 1;
    }

    const int64_t ppd64 = zeldovich_params_get_ppd(params);
    int narray = 0;
    const int narray_rc = wisdom_narray_from_params(params, &narray);
    zeldovich_params_destroy(params);
    if (narray_rc != 0) {
        MPI_Barrier(MPI_COMM_WORLD);
        return 1;
    }

    if (ppd64 <= 0 || ppd64 > (int64_t)INT_MAX) {
        fprintf(stderr, "IC_Rank0Wisdom: invalid ppd=%lld from %s\n", (long long)ppd64, param_file);
        MPI_Barrier(MPI_COMM_WORLD);
        return 1;
    }

    const int N = (int)ppd64;
    const size_t nbytes = (size_t)N * (size_t)N * sizeof(fftw_complex_t);
    fftw_complex_t *plan_buffer = NULL;
    if (posix_memalign((void **)&plan_buffer, ALIGN_BYTES, nbytes) != 0) {
        fprintf(stderr, "IC_Rank0Wisdom: posix_memalign failed (%zu bytes)\n", nbytes);
        MPI_Barrier(MPI_COMM_WORLD);
        return 1;
    }

    fftw_plan_t plan_2d = NULL;
    fftw_plan_t plan_1d = NULL;
    const int rc = wisdom_rank0_plans(N, narray, plan_buffer, &plan_2d, &plan_1d, 1);
    if (plan_2d) {
        FFTW_DESTROY_PLAN(plan_2d);
    }
    if (plan_1d) {
        FFTW_DESTROY_PLAN(plan_1d);
    }
    free(plan_buffer);

    MPI_Barrier(MPI_COMM_WORLD);
    return rc == 0 ? 0 : 1;
}

int IC_ParamBuffer(const char *bytes, size_t len, const char *param_path, const char *wisdom_save_dir)
{
    if (bytes == NULL || len < 2 || param_path == NULL || param_path[0] == '\0') {
        fprintf(stderr, "IC_ParamBuffer: invalid header bytes or param_path\n");
        return 1;
    }

    zd_ic_wisdom_save_dir = wisdom_save_dir;

    const int wis_rc = wisdom_preflight_from_param_buffer(bytes, len, param_path, wisdom_save_dir);
    if (wis_rc != 0) {
        return wis_rc;
    }

    char prog[] = "Zeldovich_MPI";
    std::string path_copy(param_path);
    char *argv[] = {prog, path_copy.data(), NULL};

    zeldovich_embed_param_header.bytes = bytes;
    zeldovich_embed_param_header.len = len;
    const int driver_rc = zeldovich_mpi_driver_run(2, argv);
    zeldovich_embed_param_header.bytes = NULL;
    zeldovich_embed_param_header.len = 0;
    return driver_rc;
}

int IC_Run(int argc, char **argv)
{
    return zeldovich_mpi_driver_run(argc, argv);
}

void IC_FinalizeStage(void)
{
    MPI_Barrier(MPI_COMM_WORLD);
    zeldovich_ic_embedded = false;
    zeldovich_embed_param_header.bytes = NULL;
    zeldovich_embed_param_header.len = 0;
    zd_ic_wisdom_save_dir = NULL;
    zd_wisdom_set_preflight_broadcast_done(false);
    zd_wisdom_set_dir(NULL);
}

} // extern "C"

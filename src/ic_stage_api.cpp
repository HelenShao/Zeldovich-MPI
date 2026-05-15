#include "ic_stage_api.h"

#include "config.h"
#include <limits.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "fft/wisdom_rank0.h"
#include "ic_embed_flags.h"
#include "zeldovich_mpi_driver.h"
#include "zeldovich_wrapper.h"

bool zeldovich_ic_embedded = false;

extern "C" {

void IC_InitStage(int from_abacus_host)
{
    // so the driver knows Abacus owns MPI / embed rules apply.
    zeldovich_ic_embedded = (from_abacus_host != 0);
}

int IC_Rank0Wisdom(const char *param_file)
{
    /* 
    Get N / narray from the parameter file, 
    allocate the buffer, call wisdom_rank0_plans_and_export,
    then tear down plans and free the buffer.
    */
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
    const int qdensity = zeldovich_params_get_qdensity(params);
    int narray;
    if (qdensity == 2) {
        narray = 1;
    } else {
        const int qPLT = zeldovich_params_get_qPLT(params);
        narray = qPLT ? 4 : 2;
    }
    zeldovich_params_destroy(params);

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
    const int rc = wisdom_rank0_plans_and_export(N, narray, plan_buffer, &plan_2d, &plan_1d);
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

int IC_Run(int argc, char **argv)
{
    // same CLI entry as standalone, but MPI already initialized by Abacus.
    return zeldovich_mpi_driver_run(argc, argv);
}

void IC_FinalizeStage(void)
{
    MPI_Barrier(MPI_COMM_WORLD); // sync all ranks after IC
    zeldovich_ic_embedded = false; // so later code does not think it is still in the embedded IC phase.
}

} // extern "C"

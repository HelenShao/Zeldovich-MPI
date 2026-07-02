/* This module owns creation, querying, and teardown of zd_comm_2d*/

#include "mpi_topology.h"

#include <stdio.h>
#include <stdlib.h>

#include "ic_embed_flags.h"

MPI_Comm zd_comm_2d = MPI_COMM_NULL; // 2D ZD cart
int zd_cart_size_x = 0; // Cart dimensions — Abacus naming: dim 0 = x-split count, dim 1 = NumZRanks
int zd_cart_size_z = 0;
int zd_cart_rank_x = 0; // coords[0] and coords[1] from the cart
int zd_cart_rank_z = 0;
int zd_cart_rank = 0; // This process’s rank in zd_comm_2d (may differ from MPI_COMM_WORLD rank when reorder=1)
int zd_comm_2d_owned = 0; // If true, zd_topology_cleanup() may MPI_Comm_free the dup/created comm

static int zd_topology_query_cart(void)
{
    // Private helper run after the comm exists. 
    // Updates zd_cart_size_x, zd_cart_size_z, zd_cart_rank, zd_cart_rank_x, zd_cart_rank_z
    if (zd_comm_2d == MPI_COMM_NULL) {
        return 1;
    }

    int ndims = 0;
    int dims[2] = {0, 0};
    int periodic[2] = {0, 0};
    int coords_unused[2] = {0, 0};
    MPI_Cart_get(zd_comm_2d, 2, dims, periodic, coords_unused);
    zd_cart_size_x = dims[0];
    zd_cart_size_z = dims[1];

    MPI_Comm_rank(zd_comm_2d, &zd_cart_rank);

    int coords[2] = {0, 0};
    MPI_Cart_coords(zd_comm_2d, zd_cart_rank, 2, coords);
    zd_cart_rank_x = coords[0];
    zd_cart_rank_z = coords[1];

    return 0;
}

int zd_topology_init_embedded(MPI_Comm abacus_comm_2d, int expected_size_z)
{
    // Called from zeldovich_mpi_driver.cpp when zeldovich_ic_embedded is true.
    int world_rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (abacus_comm_2d == MPI_COMM_NULL) {
        if (world_rank == 0) {
            fprintf(stderr, "zd_topology_init_embedded: abacus_comm_2d is MPI_COMM_NULL\n");
        }
        return 1;
    }

    zd_topology_cleanup(); // reset any prior state

    int dup_err = MPI_Comm_dup(abacus_comm_2d, &zd_comm_2d);
    if (dup_err != MPI_SUCCESS || zd_comm_2d == MPI_COMM_NULL) {
        if (world_rank == 0) {
            fprintf(stderr, "zd_topology_init_embedded: MPI_Comm_dup failed (err=%d)\n", dup_err);
        }
        return 1;
    }
    zd_comm_2d_owned = 1;
    MPI_Comm_set_errhandler(zd_comm_2d, MPI_ERRORS_RETURN);

    // cart query - if 1 then it failed, so we clean up and return error
    if (zd_topology_query_cart() != 0) {
        zd_topology_cleanup();
        return 1;
    }

    int num_ranks = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

    if (zd_cart_size_z != expected_size_z) {
        if (world_rank == 0) {
            fprintf(stderr,
                    "zd_topology_init_embedded: cart size_z=%d != NumZRanks=%d\n",
                    zd_cart_size_z,
                    expected_size_z);
        }
        zd_topology_cleanup();
        return 1;
    }
    if (zd_cart_size_x * zd_cart_size_z != num_ranks) {
        if (world_rank == 0) {
            fprintf(stderr,
                    "zd_topology_init_embedded: size_x*size_z=%d*%d != num_ranks=%d\n",
                    zd_cart_size_x,
                    zd_cart_size_z,
                    num_ranks);
        }
        zd_topology_cleanup();
        return 1;
    }

    if (zd_cart_rank == 0) {
        printf("========================================================================\n");
        printf("MPI Cartesian Topology (embedded Abacus comm_2d dup)\n");
        printf("  Grid: size_x=%d x size_z=%d = %d ranks (NumZRanks=size_z=%d)\n",
               zd_cart_size_x,
               zd_cart_size_z,
               num_ranks,
               expected_size_z);
        printf("  rank_x=coords[0], rank_z=coords[1] (Abacus-aligned)\n");
        printf("========================================================================\n");
    }

    printf("[Rank %d] Cartesian coords: (rank_x=%d, rank_z=%d)\n",
           zd_cart_rank,
           zd_cart_rank_x,
           zd_cart_rank_z);

    return 0;
}

int zd_topology_init_standalone(int num_ranks, int num_z_ranks)
{
    int world_rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (num_z_ranks <= 0 || num_ranks % num_z_ranks != 0) {
        if (world_rank == 0) {
            fprintf(stderr,
                    "zd_topology_init_standalone: num_ranks=%d not divisible by NumZRanks=%d\n",
                    num_ranks,
                    num_z_ranks);
        }
        return 1;
    }

    zd_topology_cleanup();

    // was: grid_x = NumZRanks, grid_z = num_ranks / grid_x; dims = {grid_z, grid_x}
    const int size_z = num_z_ranks;
    const int size_x = num_ranks / size_z;

    // --- Phase 7 (Abacus-aligned naming) ---
    int dims[2] = {size_x, size_z};
    int periodic[2] = {1, 1};
    const int reorder = 1;

    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periodic, reorder, &zd_comm_2d);
    if (zd_comm_2d == MPI_COMM_NULL) {
        if (world_rank == 0) {
            fprintf(stderr,
                    "zd_topology_init_standalone: MPI_Cart_create failed size_x=%d size_z=%d\n",
                    size_x,
                    size_z);
        }
        return 1;
    }
    zd_comm_2d_owned = 1;
    MPI_Comm_set_errhandler(zd_comm_2d, MPI_ERRORS_RETURN);

    if (zd_topology_query_cart() != 0) {
        zd_topology_cleanup();
        return 1;
    }

    if (zd_cart_rank == 0) {
        printf("========================================================================\n");
        printf("MPI Cartesian Topology Initialized (standalone, parameter-driven)\n");
        printf("  Grid: size_x=%d x size_z=%d = %d ranks (NumZRanks=size_z=%d)\n",
               zd_cart_size_x,
               zd_cart_size_z,
               num_ranks,
               num_z_ranks);
        printf("  Periodic: [X=%s, Z=%s]\n",
               periodic[0] ? "yes" : "no",
               periodic[1] ? "yes" : "no");
        printf("  Reorder: %s (hardware-aware rank assignment)\n", reorder ? "enabled" : "disabled");
        printf("========================================================================\n");
    }

    printf("[Rank %d] Cartesian coords: (rank_x=%d, rank_z=%d)\n",
           zd_cart_rank,
           zd_cart_rank_x,
           zd_cart_rank_z);

    return 0;
}

void zd_topology_cleanup(void)
{
    if (zd_comm_2d_owned && zd_comm_2d != MPI_COMM_NULL) {
        MPI_Comm_free(&zd_comm_2d);
    }
    zd_comm_2d = MPI_COMM_NULL;
    zd_comm_2d_owned = 0;
    zd_cart_size_x = 0;
    zd_cart_size_z = 0;
    zd_cart_rank_x = 0;
    zd_cart_rank_z = 0;
    zd_cart_rank = 0;
    zd_abacus_host_comm_2d = MPI_COMM_NULL;
}

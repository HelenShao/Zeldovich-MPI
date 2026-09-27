#ifndef MPI_TOPOLOGY_H
#define MPI_TOPOLOGY_H

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Abacus-aligned 2D Cartesian grid (Phase 7):
 *   dims[0] = size_x  (cart dim 0 = Abacus MPI_size_x)
 *   dims[1] = size_z  (cart dim 1 = Abacus MPI_size_z = NumZRanks)
 *   rank_x = coords[0], rank_z = coords[1]
 *   linear rank = rank_x * size_z + rank_z
 *
 * was: dims [grid_x, grid_z] with transposed naming (grid_x = NumZRanks)
 *
 * Embedded: zd_comm_2d = MPI_Comm_dup(Abacus comm_2d) passed via IC_InitStage.
 * Standalone: MPI_Cart_create(MPI_COMM_WORLD, dims={size_x, size_z}, reorder=1).
 */

extern MPI_Comm zd_comm_2d;
extern int zd_cart_size_x;
extern int zd_cart_size_z;
extern int zd_cart_rank_x;
extern int zd_cart_rank_z;
extern int zd_cart_rank;

/** Non-zero when Zeldovich must MPI_Comm_free zd_comm_2d in zd_topology_cleanup. */
extern int zd_comm_2d_owned;

/** Dup Abacus comm_2d; validate size_z == expected_size_z (NumZRanks). */
int zd_topology_init_embedded(MPI_Comm abacus_comm_2d, int expected_size_z);

/** Create cart on MPI_COMM_WORLD; size_z = NumZRanks, size_x = num_ranks / size_z. */
int zd_topology_init_standalone(int num_ranks, int num_z_ranks);

/** Free zd_comm_2d only when owned; reset cart globals. */
void zd_topology_cleanup(void);

#ifdef __cplusplus
}
#endif

#endif

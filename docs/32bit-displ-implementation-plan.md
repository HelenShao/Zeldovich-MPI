# Implementation Plan: Replace MPI_Ialltoallv with MPI_Isend/MPI_Irecv

## Purpose

Fix the 32-bit displacement overflow in `MPI_Ialltoallv` when `recv_total_elems ≥ 2³¹` (e.g., N=4096, narray=4, num_ranks=128). Also future-proof for N=32,000 and 81×81 nodes.

---

## Clarifications (Questions Answered)

### 1. What is the purpose of `MPI_Request *requests`? Is it to keep track of the number of sends/recvs?

`MPI_Request` is not a count—it is an MPI handle for each non-blocking operation. When you call `MPI_Isend` or `MPI_Irecv`, MPI returns an `MPI_Request` for that specific send/recv. You must pass all these handles to `MPI_Waitall` so MPI knows which operations to wait on. So `requests[]` holds the *handles* for every posted Isend and Irecv; `nreqs` is the count. Both are needed.

### 2. What is zero-count peers?

A *peer* is another rank. *Zero-count* means `sendcounts_batch[dest] == 0` or `recvcounts_batch[src] == 0` for that peer. For example, in some batches a rank may have no data to send to certain destinations (e.g., idle ranks, or when a rank has no slices in that batch). Posting Isend/Irecv with count 0 is unnecessary and adds overhead. We only post for peers where the count is > 0.

### 3. Receives first: Registers receive buffers before any sends start, avoiding internal buffering pressure on the MPI library

When you call `MPI_Isend`, MPI may need to buffer the outgoing data until the matching `MPI_Irecv` is posted. If you post all sends first, the MPI library must hold all that data internally, which can exhaust buffers and cause blocking or failures. By posting all `Irecv` first, every send has a ready receive buffer, so data flows directly send→recv without heavy internal buffering. This order is a common best practice.

### 4. Do I still need an MPI Cartesian grid for the MPI decomposition?

Yes. The 2D Cartesian topology (`comm_2d`) is used for the XZ pencil decomposition: it defines which (X,Z) region each rank owns. The Isend/Irecv change only replaces *how* the all-to-all data exchange is done, not *which* data is sent or where it goes. The decomposition logic, `get_padded_bounds_simple`, and all grid layout remain unchanged.

### 5. With Irecv/Isend, how do I unpack from the recv buffer?

The same way as with Alltoallv. The recv buffer layout is identical: data from source rank `src` for this batch is written at `recv_buffer + rdispls_elem[src]`, and the size is `recvcounts_batch[src]` elements. The streaming stage (`z_streaming_unpack`) uses `recv_displs_src` and `src_write_cursor` to compute offsets; those are unchanged. Isend/Irecv only changes *how* MPI fills the buffer, not *where* it puts the data.

### 6. Do `rdispls_elem` and `sdispls_elem` avoid overflow because they use 64-bit int?

Yes. `int64_t` can represent values up to 2⁶³ − 1. `rdispls_elem[src]` and `sdispls_elem[dest]` are element offsets that can exceed `INT_MAX` (2³¹ − 1) when the total buffer is large. By using `int64_t`, we avoid overflow. We then pass `recv_buffer + rdispls_elem[src]` (a 64-bit pointer) to `MPI_Irecv`, so MPI never sees a 32-bit displacement—only a pointer. Same for sends.

### 7. Modify `calculate_batch_send_recv_counts` to return `int64_t` displacements and totals

Included in the plan below. The function will return `int64_t *sdispls`, `int64_t *rdispls`, and `int64_t *total_send`, `int64_t *total_recv` to support N=32,000.

### 8. Widen `pack_slices_to_send_buffer` internal indexing to `int64_t`

Included in the plan below. The function will accept `int64_t *sdispls` and use `int64_t` for `total_send_size`, `offset`, and buffer indexing.

---

## Overview

The change is localized to the batch communication step (Steps 3–7 in the batch loop) plus updates to `batch_helpers.c` and `pack_slices_to_send_buffer`. Everything else (generation, streaming, unpacking) stays the same.

---

## File 1: `src/main.cpp` — Batch loop (lines ~605–790)

### Step 3: Use int64_t displacements, remove INT_MAX overflow checks

- Remove the `#if 1` block that checks for INT_MAX overflow and fills `rdispls_batch[src]`
- Keep `rdispls_elem[src]` = `recv_displs_src[src] + src_write_cursor[src]` (already int64_t)
- Compute `sdispls_elem[dest]` from `calculate_batch_send_recv_counts` (it will now return int64_t sdispls)
- Remove the two `#if 1` blocks that abort when displacements exceed INT_MAX
- Keep the bounds-checking logic (rdispls_elem + recvcounts vs recv_buffer_size)

### Step 6: Replace MPI_Ialltoallv with Isend/Irecv loop

Replace:

```cpp
MPI_Request comm_request_batch;
MPI_Ialltoallv(send_buffer_batch, sendcounts_batch, sdispls_batch, MPI_COMPLEX_TYPE,
               recv_buffer, recvcounts_batch, rdispls_batch, MPI_COMPLEX_TYPE,
               comm_2d, &comm_request_batch);
MPI_Wait(&comm_request_batch, MPI_STATUS_IGNORE);
```

With:

```cpp
MPI_Request *requests = (MPI_Request *)malloc(sizeof(MPI_Request) * 2 * num_ranks);
int nreqs = 0;

// Post all receives first (recommended ordering to avoid MPI internal buffering)
for (int src = 0; src < num_ranks; src++) {
    if (recvcounts_batch[src] > 0) {
        MPI_Irecv(recv_buffer + rdispls_elem[src],
                  recvcounts_batch[src], MPI_COMPLEX_TYPE,
                  src, batch_idx, comm_2d, &requests[nreqs++]);
    }
}

// Then post all sends
for (int dest = 0; dest < num_ranks; dest++) {
    if (sendcounts_batch[dest] > 0) {
        MPI_Isend(send_buffer_batch + sdispls_elem[dest],
                  sendcounts_batch[dest], MPI_COMPLEX_TYPE,
                  dest, batch_idx, comm_2d, &requests[nreqs++]);
    }
}

MPI_Waitall(nreqs, requests, MPI_STATUSES_IGNORE);
free(requests);
```

Design decisions:
- **Tag = `batch_idx`**: Ensures messages from different batches do not match
- **Receives first**: Registers receive buffers before sends start
- **Zero-count peers**: Skip Isend/Irecv when count is 0
- **Self-send**: MPI handles rank sending to itself when both are non-blocking

### Step 7

MPI_Wait is folded into MPI_Waitall above; remove the separate `MPI_Wait`.

### Step 9: Cleanup

After modifying `calculate_batch_send_recv_counts`, it will return `int64_t *sdispls` and `int64_t *rdispls`. Update frees accordingly. `rdispls_batch` (from the old API) is no longer needed.

---

## File 2: `src/utils/batch_helpers.c` — Return int64_t displacements and totals

### Update `calculate_batch_send_recv_counts` signature

Change from:

```c
void calculate_batch_send_recv_counts(
    int rank, int num_ranks, int N, int narray, int batch_idx,
    int my_batch_slice_count, int my_pencils,
    int **out_sendcounts, int **out_sdispls,
    int **out_recvcounts, int **out_rdispls,
    int *out_total_send, int *out_total_recv);
```

To:

```c
void calculate_batch_send_recv_counts(
    int rank, int num_ranks, int N, int narray, int batch_idx,
    int my_batch_slice_count, int my_pencils,
    int **out_sendcounts, int64_t **out_sdispls,
    int **out_recvcounts, int64_t **out_rdispls,
    int64_t *out_total_send, int64_t *out_total_recv);
```

### Implementation changes

- Allocate `int64_t *sdispls` and `int64_t *rdispls` instead of `int *`
- Use `int64_t total_send` and `int64_t total_recv` in the accumulation loops
- `sendcounts` and `recvcounts` stay `int` (per-peer counts fit in int even at N=32,000)

### Update `batch_helpers.h` accordingly

---

## File 3: `src/communication/mpi_exchange.c` — Widen pack_slices_to_send_buffer

### Update `pack_slices_to_send_buffer` signature

Change from:

```c
void pack_slices_to_send_buffer(
    int rank, int num_ranks, int N, int narray,
    fftw_complex_t *local_y_slices,
    int num_my_slices, int *y_global_map,
    fftw_complex_t *send_buffer, int *sendcounts, int *sdispls);
```

To:

```c
void pack_slices_to_send_buffer(
    int rank, int num_ranks, int N, int narray,
    fftw_complex_t *local_y_slices,
    int num_my_slices, int *y_global_map,
    fftw_complex_t *send_buffer, int *sendcounts, int64_t *sdispls);
```

### Implementation changes

- `int offset` → `int64_t offset`
- `int total_send_size` → `int64_t total_send_size`
- `int dest_offset` → `int64_t dest_offset`
- All buffer index calculations use `int64_t` to avoid overflow at N=32,000
- Update `mpi_exchange.h` accordingly

---

## File 4: `src/config.h` — Optional config flag

```c
// Communication method for all-to-all exchange
// 0 = MPI_Ialltoallv (requires 32-bit displacements, may overflow)
// 1 = MPI_Isend/MPI_Irecv (no displacement limit, recommended)
#ifndef USE_ISEND_IRECV
#define USE_ISEND_IRECV 1
#endif
```

Use `#if USE_ISEND_IRECV ... #else ... #endif` around the communication block if both paths are desired for comparison. Otherwise, fully switch to Isend/Irecv.

---

## Summary of Changes

| Component | Changes |
|-----------|---------|
| `main.cpp` | Replace Alltoallv with Isend/Irecv; use int64_t sdispls_elem, rdispls_elem; remove INT_MAX guards |
| `batch_helpers.c` | Return int64_t sdispls, rdispls, total_send, total_recv |
| `batch_helpers.h` | Update function signature |
| `mpi_exchange.c` | Accept int64_t sdispls; use int64_t for internal indexing |
| `mpi_exchange.h` | Update pack_slices_to_send_buffer signature |
| `config.h` | Add USE_ISEND_IRECV flag (optional) |
| `z_streaming.c` | No changes (already uses int64_t recv_displs_src) |

---

## What Stays the Same

- MPI Cartesian topology and 2D grid decomposition
- Generation, packing layout, streaming, unpacking logic
- recv_buffer layout and recv_displs_src, src_write_cursor usage

---

## Testing Plan

1. **N=4, num_ranks=2**: Sanity check
2. **N=1024, num_ranks=9 or 32**: Compare output to existing verified runs
3. **N=4096, num_ranks=128**: Previously overflowed; should run successfully
4. **N=32,000, 81×81 nodes**: Target scale; verify completion and correctness
5. **Performance**: Compare Isend/Irecv vs Alltoallv at N=1024

---

## Future: Overlap (Optional)

To overlap communication of batch N with generation of batch N+1:
- Double-buffer send buffers
- Start comm for batch 0; for batch 1: generate+pack; Wait(comm 0); start comm 1; etc.
- Same pattern works with Isend/Irecv.

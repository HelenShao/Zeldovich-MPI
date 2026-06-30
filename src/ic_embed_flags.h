#pragma once

#include <stddef.h>

/** True when Zeldovich_IC runs inside an Abacus host (skip global FFTW cleanup, etc.). */
extern bool zeldovich_ic_embedded;

/** When set, zeldovich_mpi_driver_run uses these bytes instead of re-reading the param file. */
struct ZeldovichEmbedParamHeader {
    const char *bytes;
    size_t len;
};
extern ZeldovichEmbedParamHeader zeldovich_embed_param_header;

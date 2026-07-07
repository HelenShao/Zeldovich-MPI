#pragma once

#include <fmt/core.h>

#ifdef __cplusplus
extern "C" {
#endif

/** True on MPI world rank 0 (informational logs should use this in MPI runs). */
bool zd_log_on_this_rank(void);

#ifdef __cplusplus
}
#endif

#define ZD_ERR(...) \
    do { \
        if (zd_log_on_this_rank()) { \
            fmt::print(stderr, __VA_ARGS__); \
        } \
    } while (0)

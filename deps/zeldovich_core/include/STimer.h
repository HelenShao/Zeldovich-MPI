#ifndef STIMER
#define STIMER

#include <cassert>
#include <time.h>

// N.B. STimer has no cache line padding because some containers
// (such as those in TBB) use std::aligned_storage internally, which
// may not support over-alignment.

class STimer {
public:
    STimer();
    ~STimer();
    struct timespec Start(void);
    struct timespec StartFrom(struct timespec);
    struct timespec Stop(void);
    void StopAt(struct timespec);
    double Elapsed(void) const;
    void Clear(void);
    void increment(struct timespec dt);
    struct timespec get_timer(void) const;
    struct timespec get_start(void) const;

    int timeron;
    struct timespec timer;

private:
    struct timespec tstart;
};

struct timespec scale_timer(double s, struct timespec t);
void timespecclear(struct timespec *t);
void timespecadd(struct timespec *a, struct timespec *b, struct timespec *res);
void timespecsub(struct timespec *a, struct timespec *b, struct timespec *res);
#endif

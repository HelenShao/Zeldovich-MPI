#ifndef KERNEL_TERMS_H
#define KERNEL_TERMS_H

// ====================================================================================
// KERNEL TERM TABLES (f_NL kernels for Kernel_MPI)
// ====================================================================================
// Every kernel used by kernel_MPI.py has the form
//     G(k) = |k|^pG,  H(k) = |k|^pH,  J(k) = cJ * |k|^pJ
// and contributes  coef * J(k) * Q(k),  Q = FFT[ IFFT(G*phi) * IFFT(H*phi) ].
// Each f_NL type is a list of (shape, coef) terms summed into one output.
// ====================================================================================

#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    KSHAPE_A1 = 0,  // G=1,     H=1,     J=1
    KSHAPE_C1,      // G=k,     H=k,     J=1/(3k^2)
    KSHAPE_B1,      // G=k,     H=1,     J=1/k
    KSHAPE_B2,      // G=k^2,   H=1,     J=1/k^2
    KSHAPE_D1,      // G=1/k,   H=1/k,   J=k^2
    KSHAPE_E2,      // G=k^2,   H=1/k,   J=1/(2k)
    KSHAPE_F1,      // G=1/k,   H=1,     J=k
    KSHAPE_G1,      // G=1/k,   H=k,     J=1
    KSHAPE_COUNT
} KernelShapeId;

typedef struct {
    const char *name;
    int pG, pH, pJ;   // powers of |k| for G, H, J
    double cJ;        // constant prefactor of J
} KernelShape;

typedef struct {
    KernelShapeId shape;
    double coef;      // coefficient of this term in the f_NL sum
} KernelTerm;

typedef enum {
    FNL_NONE = -1,
    FNL_LOCAL = 0,
    FNL_EQUIL,
    FNL_ORTHOCMB,
    FNL_ORTHOLSS
} FnlType;

#define KERNEL_MAX_TERMS 8

const KernelShape *kernel_shape(KernelShapeId id);

// Returns the shape id for a name like "A1" or "C1", or -1 if unknown
int kernel_shape_from_name(const char *name);

const char *fnl_type_name(FnlType type);

// Fills terms[] (size >= KERNEL_MAX_TERMS) in the same order kernel_MPI.py adds them.
// Returns the number of terms.
int kernel_terms_for(FnlType type, KernelTerm *terms);

// |k|^p with the kernel_MPI.py DC convention: kpow(0, 0) = 1, kpow(0, p != 0) = 0
static inline double kpow(double k, int p)
{
    if (p == 0) return 1.0;
    if (k == 0.0) return 0.0;
    switch (p) {
        case 1:  return k;
        case 2:  return k * k;
        case -1: return 1.0 / k;
        case -2: return 1.0 / (k * k);
        default: return pow(k, (double)p);
    }
}

#ifdef __cplusplus
}
#endif

#endif

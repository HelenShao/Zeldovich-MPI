#pragma once

#include <stdint.h>

#include <vector>

#include <pcg-rng/pcg_random.hpp>

#ifdef HAVE_GSL
#include <gsl/gsl_rng.h>
#endif

#include <stddef.h>

#include "parameters.h"
#include "spline_function.h"
#include "zeldovich.h"

/// Rank 0: load P(k) file into vectors.
int ReadPkTextFileIntoVectors(const fs::path &filename, Parameters &param,
    std::vector<double> &k_out, std::vector<double> &p_out);

class PowerSpectrum : public SplineFunction {
public:
    PowerSpectrum(int n, Parameters &param);
    ~PowerSpectrum();

    int fixed_power;
    int is_powerlaw;
    double powerlaw_index;
    double normalization;
    double Pk_smooth2;  // param.Pk_smooth squared
    double Rnorm;
    double kmax;  // max k in the input PS
    double kmin;  // min (non-zero) k in the input PS
    int block;         // ppd/numblock — v1 slab height only; NOT v2rng length
    int v2rng_count;   // ppd/2 when v2rng is active (global Y index range)
    double primordial_norm;
    double n_s;

#ifdef HAVE_GSL
    gsl_rng **v1rng;  // The random number generators for the deprecated ZD_Version=1
#endif
    pcg64 *v2rng;     // The random number generators for version 2 (current version)

    double sigmaR_integrand(double k);
    double sigmaR(double R);

    double Romberg(
       double (PowerSpectrum::*func)(double),
       double a,
       double b,
       double prec,
       double *obtprec
    );

    int InitFromFile(const fs::path &filename, Parameters &param);

    // builds spline-interpolated pk from arrays that are in memory, instead of reading from file
    int InitFromRawPk(const double *k_arr, const double *p_arr, size_t n, Parameters &param);

    int InitFromPowerLaw(double _powerlaw_index, Parameters &param);

    void Normalize(Parameters &param);

    double power(double wavenumber);
    double primordial_power(double wavenumber);
    double infer_Tk(double wavenumber);

    template <int Ver>
    double one_rand(int64_t i);

    template <int Ver>
    Complx cgauss(double wavenumber, int64_t rng);
    
    // Get a copy of the RNG for specific Y-slice (for thread-local use)
    pcg64 get_rng_copy(int64_t rng_index) const;
};

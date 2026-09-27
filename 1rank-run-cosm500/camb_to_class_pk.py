"""
Compute the late-time (z=0) and primordial matter power spectra with CLASS,
using the same cosmology as the CAMB run in camb_planck15_params.ini, and
write them out in the same two-column [k, P(k)] format as camb_planck15.dat.

Output files:
    class_pk_z0.dat         -> k [h/Mpc],  P(k,z=0) [(Mpc/h)^3]   (late-time)
    class_pk_primordial.dat -> k [h/Mpc],  P_R(k)   [dimensionless] (primordial)

Note on units/interpretation:
    The late-time file is directly comparable to camb_planck15.dat (matterpower.dat):
    same units, same convention (k in h/Mpc, P(k) in (Mpc/h)^3).

    The primordial file is NOT in the same physical units as the matter power
    spectrum -- P_R(k) is the dimensionless power spectrum of the curvature
    perturbation, P_R(k) = A_s (k/k_pivot)^(n_s-1), before any transfer function
    is applied. It's written with the same k-grid (h/Mpc) purely so the two
    files line up column-for-column for easy plotting/comparison; the second
    column is a different physical quantity, not a matter power spectrum.
"""

import numpy as np
from classy import Class

# ---------------------------------------------------------------------
# 1) Read the k-grid from the CAMB output, so both CLASS files use the
#    exact same k values as camb_planck15.dat for direct comparison.
# ---------------------------------------------------------------------
camb_data = np.loadtxt("data/camb_planck15.dat")
k_hmpc = camb_data[:, 0]          # k in h/Mpc
kmin_hmpc, kmax_hmpc = k_hmpc.min(), k_hmpc.max()

# ---------------------------------------------------------------------
# 2) Map the CAMB ini parameters onto the equivalent CLASS parameters.
#
#    ombh2        = 0.02222   -> omega_b
#    omch2        = 0.1199    -> omega_cdm
#    omnuh2       = 0.0, massive_neutrinos = 0 -> no massive neutrinos
#    massless_neutrinos = 3.04 -> N_ur = 3.04 (CLASS's massless-neutrino species)
#    hubble       = 67.26     -> h = 0.6726
#    temp_cmb     = 2.7255    -> T_cmb
#    helium_fraction = 0.2453 -> YHe
#    scalar_amp(1)        = 2.1e-9  -> A_s
#    scalar_spectral_index(1) = 0.9652 -> n_s
#    scalar_nrun(1)       = 0  -> alpha_s (no running)
#    pivot_scalar = 0.05       -> k_pivot
#    do_nonlinear = 0          -> non linear = none
# ---------------------------------------------------------------------
params = {
    'output': 'mPk',
    'non linear': 'none',
    'gauge': 'Newtonian', # compute perturbations in Newtonian (longitudinal) gauge

    'h': 0.6726,
    'omega_b': 0.02222,
    'omega_cdm': 0.1199,
    'N_ur': 3.04,
    'N_ncdm': 0,

    'T_cmb': 2.7255,
    'YHe': 0.2453,

    'A_s': 2.1e-9,
    'n_s': 0.9652,
    'alpha_s': 0.0,
    'k_pivot': 0.05,

    # Make sure CLASS tabulates P(k) over (at least) the full CAMB k-range
    'P_k_max_h/Mpc': kmax_hmpc * 1.05,
    'z_pk': 0,
}

cosmo = Class()
cosmo.set(params)
cosmo.compute()

h = cosmo.h()

# ---------------------------------------------------------------------
# 3) Late-time matter power spectrum, z = 0
#    classy's cosmo.pk(k, z) expects k in 1/Mpc and returns P(k) in Mpc^3.
#    Convert to the CAMB convention: k in h/Mpc, P(k) in (Mpc/h)^3.
# ---------------------------------------------------------------------
k_1mpc = k_hmpc * h
pk_z0_mpc3 = np.array([cosmo.pk(k, 0.0) for k in k_1mpc])   # Mpc^3
pk_z0_hmpc3 = pk_z0_mpc3 * h**3                              # (Mpc/h)^3

late_time_out = np.column_stack([k_hmpc, pk_z0_hmpc3])
np.savetxt(
    "data/class_pk_z0.dat",
    late_time_out,
    fmt="%15.8E",
)

# ---------------------------------------------------------------------
# 4) Primordial (curvature) power spectrum P_R(k)
#    get_primordial() returns CLASS's internal k-grid (1/Mpc) and P_scalar(k).
#    Interpolate (in log-log space, since it's a smooth power law) onto the
#    same k-grid as above so the two files can be compared column-by-column.
# ---------------------------------------------------------------------
prim = cosmo.get_primordial()
k_prim_1mpc = np.array(prim['k [1/Mpc]'])
p_prim = np.array(prim['P_scalar(k)'])

log_p_interp = np.interp(
    np.log(k_1mpc),
    np.log(k_prim_1mpc),
    np.log(p_prim),
)
p_prim_at_k = np.exp(log_p_interp)

# primordial_out = np.column_stack([k_hmpc, p_prim_at_k])
# np.savetxt(
#     "/mnt/user-data/outputs/class_pk_primordial.dat",
#     primordial_out,
#     fmt="%15.8E",
# )

# ---------------------------------------------------------------------
# 5) Convert the dimensionless primordial curvature spectrum P_R(k)
#    (what get_primordial() returns, = A_s (k/k*)^(n_s-1)) into a
#    volume-dimensioned power spectrum, in the SAME units as the
#    late-time file: (Mpc/h)^3, on the same k_h/Mpc grid.
#
#        P_R(k) [Mpc^3]        = (2 pi^2 / k^3) * P_scalar(k)     (k in 1/Mpc)
#        P_R(k_h) [(Mpc/h)^3]  = (2 pi^2 / k_h^3) * P_scalar(k_h*h)
#
#    (the h^3 factors cancel between the Mpc^3 -> (Mpc/h)^3 conversion
#    and the k -> k_h conversion, provided P_scalar is evaluated at the
#    *physical* k = k_h*h, since k_pivot is defined in 1/Mpc)
# ---------------------------------------------------------------------
p_prim_volume_hmpc3 = (2.0 * np.pi**2 / k_hmpc**3) * p_prim_at_k

primordial_dimensional_out = np.column_stack([k_hmpc, p_prim_volume_hmpc3])
np.savetxt(
    "data/class_pk_primordial_dimensional.dat",
    primordial_dimensional_out,
    fmt="%15.8E",
)

# Example: reconstruct delta_prim(k) from a delta_late(k) realization
# delta_prim(k) = sqrt(P_prim(k)/P_late(k,z=0)) * delta_late(k)
# (shown here just as the rescaling factor sqrt(P_prim/P_late), since
#  delta_late(k) itself is whatever complex field/realization you have)
# rescale_factor = np.sqrt(p_prim_volume_hmpc3 / pk_z0_hmpc3)
# np.savetxt(
#     "/mnt/user-data/outputs/class_pk_prim_over_late_rescale_factor.dat",
#     np.column_stack([k_hmpc, rescale_factor]),
#     fmt="%15.8E",
# )

cosmo.struct_cleanup()
cosmo.empty()

print("h =", h)
print("Wrote class_pk_z0.dat and class_pk_primordial.dat")
print("k range:", kmin_hmpc, "to", kmax_hmpc, "h/Mpc")

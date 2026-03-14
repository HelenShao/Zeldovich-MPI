#!/usr/bin/env python3
"""
Quickly print and compare the first N particle IC values from ic_0 (or any ic_* file)
for both the MPI hermitian code and zeldovich-PLT outputs.
Usage:
  python3 print_first_10_ic.py [num_particles]
  python3 print_first_10_ic.py 10
  python3 print_first_10_ic.py 10 particle_ics/ic_0 /path/to/zeldovich/ic_0 1024
  python3 print_first_10_ic.py --hermitian-only   # displacement vs velocity from hermitian only
"""

import struct
import sys
import os

BYTES_PER_PARTICLE = 32
PATTERN = '=HHH xx fff fff'  # i, j, k, padding, dx, dy, dz, vx, vy, vz
MIN_DISPL = 1e-10  # avoid division by zero for ratio vel/displ

def read_first_n_particles(filename, n, label=""):
    """Read first n particles from IC file. Returns list of (i,j,k, dx,dy,dz, vx,vy,vz) or (None, error)."""
    try:
        with open(filename, 'rb') as f:
            out = []
            for idx in range(n):
                raw = f.read(BYTES_PER_PARTICLE)
                if len(raw) < BYTES_PER_PARTICLE:
                    return out, None  # EOF, return what we got
                i, j, k, dx, dy, dz, vx, vy, vz = struct.unpack(PATTERN, raw)
                out.append((int(i), int(j), int(k), dx, dy, dz, vx, vy, vz))
            return out, None
    except FileNotFoundError:
        return None, f"File not found: {filename}"
    except Exception as e:
        return None, str(e)

def run_hermitian_only(hermitian_file, num_particles):
    """Print first N particles from hermitian ic_0 and compare displacement vs velocity."""
    herm, err = read_first_n_particles(hermitian_file, num_particles)
    if err:
        print(f"ERROR: {err}", file=sys.stderr)
        sys.exit(1)
    n_show = len(herm)
    if n_show == 0:
        print("No particles read.")
        return

    print("MPI hermitian code — first {} particles from {}".format(n_show, hermitian_file))
    print()
    print("# idx   i    j    k  |  displacement (Z, Y, X)           velocity (Z, Y, X)            |  vel/displ (Z, Y, X)  [— when |displ|<1e-10]")
    print("#" + "-" * 120)

    ratios_z, ratios_y, ratios_x = [], [], []
    for idx in range(n_show):
        (hi, hj, hk, dx, dy, dz, vx, vy, vz) = herm[idx]
        # vel order in file: vx=Z, vy=Y, vz=X (from struct: vel[0]=Z, vel[1]=Y, vel[2]=X)
        rz = vx / (dx if abs(dx) > MIN_DISPL else float('nan'))
        ry = vy / (dy if abs(dy) > MIN_DISPL else float('nan'))
        rx = vz / (dz if abs(dz) > MIN_DISPL else float('nan'))
        if abs(dx) > MIN_DISPL:
            ratios_z.append(rz)
        if abs(dy) > MIN_DISPL:
            ratios_y.append(ry)
        if abs(dz) > MIN_DISPL:
            ratios_x.append(rx)
        rz_str = "{:8.5f}".format(rz) if abs(dx) > MIN_DISPL else "      —"
        ry_str = "{:8.5f}".format(ry) if abs(dy) > MIN_DISPL else "      —"
        rx_str = "{:8.5f}".format(rx) if abs(dz) > MIN_DISPL else "      —"
        print(f"{idx:4d}  {hi:4d} {hj:4d} {hk:4d}  | "
              f"{dx:11.6f} {dy:11.6f} {dz:11.6f}   "
              f"{vx:11.6f} {vy:11.6f} {vz:11.6f}  |  {rz_str} {ry_str} {rx_str}")

    # Summary: displacement vs velocity
    print()
    print("Summary (displacement vs velocity, first {} particles):".format(n_show))
    comp = [("Z", ratios_z), ("Y", ratios_y), ("X", ratios_x)]
    for name, ratios in comp:
        if ratios:
            mean_r = sum(ratios) / len(ratios)
            print("  Component {}: vel/displ ratio  count={}, mean={:.6f}".format(name, len(ratios), mean_r))
    print("  (In linear theory with PLT, velocity ∝ f × displacement, so vel/displ ~ growth rate f.)")

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    default_hermitian = os.path.join(script_dir, "particle_ics", "ic_0")
    default_zeldovich = os.path.expanduser("~/InitialConditions/zeldovich-PLT/output_N1024_CPD_225/ic_0")
    default_n = 10

    args = [a for a in sys.argv[1:] if a in ("--hermitian-only", "-H")]
    hermitian_only = bool(args)

    if hermitian_only:
        # Consume the flag
        sys.argv = [sys.argv[0]] + [a for a in sys.argv[1:] if a not in ("--hermitian-only", "-H")]

    if hermitian_only:
        num_particles = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else default_n
        hermitian_file = sys.argv[2] if len(sys.argv) > 2 else default_hermitian
        run_hermitian_only(hermitian_file, num_particles)
        return

    if len(sys.argv) >= 5:
        num_particles = int(sys.argv[1])
        hermitian_file = sys.argv[2]
        zeldovich_file = sys.argv[3]
        N = int(sys.argv[4])
    elif len(sys.argv) == 2 and sys.argv[1].isdigit():
        num_particles = int(sys.argv[1])
        hermitian_file = default_hermitian
        zeldovich_file = default_zeldovich
        N = 1024
    else:
        num_particles = default_n
        hermitian_file = default_hermitian
        zeldovich_file = default_zeldovich
        N = 1024

    print(f"First {num_particles} particles from ic_0 (N={N})")
    print(f"  Hermitian: {hermitian_file}")
    print(f"  Zeldovich: {zeldovich_file}")
    print()

    herm, err_h = read_first_n_particles(hermitian_file, num_particles)
    if err_h:
        print(f"ERROR (hermitian): {err_h}", file=sys.stderr)
        sys.exit(1)
    zel, err_z = read_first_n_particles(zeldovich_file, num_particles)
    if err_z:
        print(f"ERROR (zeldovich): {err_z}", file=sys.stderr)
        sys.exit(1)

    n_show = min(len(herm), len(zel), num_particles)
    tol = 1e-5   # single-precision typical; use 1e-6 for strict byte-level match

    # Header
    print("# idx   i    j    k  |  Hermitian displ(Z,Y,X)              vel(Z,Y,X)              |  Zeldovich displ(Z,Y,X)              vel(Z,Y,X)              |  match")
    print("#" + "-" * 140)

    for idx in range(n_show):
        (hi, hj, hk, hdx, hdy, hdz, hvx, hvy, hvz) = herm[idx]
        (zi, zj, zk, zdx, zdy, zdz, zvx, zvy, zvz) = zel[idx]
        match_d = abs(hdx - zdx) <= tol and abs(hdy - zdy) <= tol and abs(hdz - zdz) <= tol
        match_v = abs(hvx - zvx) <= tol and abs(hvy - zvy) <= tol and abs(hvz - zvz) <= tol
        match = "OK" if (match_d and match_v) else "DIFF"
        print(f"{idx:4d}  {hi:4d} {hj:4d} {hk:4d}  | "
              f"{hdx:11.6f} {hdy:11.6f} {hdz:11.6f}  {hvx:11.6f} {hvy:11.6f} {hvz:11.6f}  | "
              f"{zdx:11.6f} {zdy:11.6f} {zdz:11.6f}  {zvx:11.6f} {zvy:11.6f} {zvz:11.6f}  | {match}")

    # Summary differences for first 10
    if n_show > 0:
        print()
        print("Summary (first {} particles):".format(n_show))
        max_dd, max_dv = 0.0, 0.0
        for idx in range(n_show):
            (hi, hj, hk, hdx, hdy, hdz, hvx, hvy, hvz) = herm[idx]
            (zi, zj, zk, zdx, zdy, zdz, zvx, zvy, zvz) = zel[idx]
            for a, b in [(hdx, zdx), (hdy, zdy), (hdz, zdz)]:
                max_dd = max(max_dd, abs(a - b))
            for a, b in [(hvx, zvx), (hvy, zvy), (hvz, zvz)]:
                max_dv = max(max_dv, abs(a - b))
        print(f"  Max |displacement difference|: {max_dd:.6e}")
        print(f"  Max |velocity difference|:     {max_dv:.6e}")
        if max_dd <= tol and max_dv <= tol:
            print("  => All values match within tolerance {:.0e}".format(tol))
        else:
            print("  => Some values differ (tolerance {:.0e})".format(tol))

if __name__ == '__main__':
    main()

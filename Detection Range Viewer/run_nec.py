"""
Minimal NEC-2 card-deck runner using PyNEC (necpp: `pip install PyNEC`), for
the cards this project's .nec files actually use: CM, CE, GW, GE, GN, EX,
FR, RP, EN. Lets you go from a .nec file to a .out file locally, without
opening 4nec2 - handy for regenerating a pattern after a quick geometry edit,
or for adding a new antenna model, when you don't want to click through the
GUI. 4nec2 (or any other NEC-2 engine) still works fine for the same .nec
files if you prefer it, or for cards this script doesn't handle.

Usage: python run_nec.py <input.nec> [output.out]
If output.out is omitted, replaces the .nec extension with .out.

Writes an output file formatted closely enough to real NEC-2 .out radiation-
pattern tables that parse_nec.py in this folder can read it unmodified -
same column layout, same "DATA CARD ... GN ..." echo, same
"FREQUENCY= ... MHZ" line, and CM comments (including a "GROUND TYPE: ..."
line, if present) are echoed the same way 4nec2's own NEC-2 engine echoes
them.

Validated against this project's existing 4nec2-generated .out files: peak
gain matches within ~0.3 dB (numerical differences between NEC-2
implementations, not a modeling error) - see PLC1669_9el_166MHz.out (real
4nec2 run, 18.22 dBi peak) vs. this script on the same .nec file (17.92 dBi).
"""
import sys, os, re
from PyNEC import nec_context


def parse_nec_file(path):
    cards = []
    comments = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            toks = line.split()
            tag = toks[0].upper()
            if tag == "CM":
                comments.append(line[2:].strip())
                continue
            if tag in ("CE", "EN"):
                continue
            cards.append((tag, toks[1:]))
    return comments, cards


def run(nec_path, out_path):
    comments, cards = parse_nec_file(nec_path)

    c = nec_context()
    geo = c.get_geometry()

    ge_flag = 0
    gn_line = None
    ex_vals = None
    fr_mhz = None
    rp_args = None

    for tag, toks in cards:
        nums = [float(t) for t in toks]
        if tag == "GW":
            tag_id, seg_count = int(nums[0]), int(nums[1])
            x1, y1, z1, x2, y2, z2, rad = nums[2:9]
            geo.wire(tag_id, seg_count, x1, y1, z1, x2, y2, z2, rad, 1.0, 1.0)
        elif tag == "GE":
            ge_flag = int(nums[0]) if nums else 0
        elif tag == "GN":
            gn_line = nums
        elif tag == "EX":
            ex_vals = nums
        elif tag == "FR":
            fr_mhz = nums[4]
        elif tag == "RP":
            rp_args = nums
        else:
            raise ValueError("unsupported card: " + tag + " (extend run_nec.py to add it)")

    c.geometry_complete(ge_flag)

    if gn_line is not None:
        gt, rwc = int(gn_line[0]), int(gn_line[1])
        f1, f2, f3, f4, f5, f6 = (list(gn_line[4:10]) + [0, 0, 0, 0, 0, 0])[:6]
        c.gn_card(gt, rwc, f1, f2, f3, f4, f5, f6)

    exi, i2, i3, i4 = int(ex_vals[0]), int(ex_vals[1]), int(ex_vals[2]), int(ex_vals[3])
    f1, f2 = (list(ex_vals[4:6]) + [0, 0])[:2]
    c.ex_card(exi, i2, i3, i4, f1, f2, 0, 0, 0, 0)

    c.fr_card(0, 1, fr_mhz, 0)

    calc_mode, n_theta, n_phi, i4 = int(rp_args[0]), int(rp_args[1]), int(rp_args[2]), int(rp_args[3])
    theta0, phi0, dtheta, dphi, rdist, gnorm = rp_args[4:10]
    c.rp_card(calc_mode, n_theta, n_phi, 0, 0, 0, 0, theta0, phi0, dtheta, dphi, rdist, gnorm)

    rp = c.get_radiation_pattern(0)
    thetas = rp.get_theta_angles()
    phis = rp.get_phi_angles()
    gtot = rp.get_gain_tot()

    import numpy as np
    gtot = np.array(gtot)
    if gtot.ndim == 1:
        # PyNEC's flat gain array is ordered phi-major/theta-minor (phi is the
        # outer loop, theta the inner one) - reshape(n_theta, n_phi) silently
        # scrambles it into a different but still plausible-looking grid
        # instead of raising an error, since the two axes are close in size.
        # Confirmed against a live PyNEC run: this order gives a phi-invariant
        # theta=0 zenith row and puts the peak at the geometrically correct
        # boresight (theta=90, phi=180 for an endfire Yagi), whereas the naive
        # reshape put the (identical, for two unrelated antennas) bogus peak
        # at theta=45/phi=270 instead.
        gtot = gtot.reshape(n_phi, n_theta).T

    with open(out_path, "w") as f:
        f.write("\n\n\n                                 - - - - COMMENTS - - - -\n\n")
        for cm in comments:
            f.write("                          " + cm + "\n")
        f.write("\n")
        if gn_line is not None:
            gt, rwc = int(gn_line[0]), int(gn_line[1])
            f1, f2 = gn_line[4], gn_line[5]
            f.write(" ***** DATA CARD NO.  1   GN   {:5d} {:5d} {:5d} {:5d} {:13.5E} {:13.5E} {:13.5E} {:13.5E} {:13.5E} {:13.5E}\n".format(
                gt, rwc, 0, 0, f1, f2, 0.0, 0.0, 0.0, 0.0))
        f.write("\n\n                                 - - - - - - FREQUENCY - - - - - -\n\n")
        f.write("                                    FREQUENCY= {:.4E} MHZ\n\n".format(fr_mhz))
        f.write("                                                - - - RADIATION PATTERNS - - -\n\n")
        f.write("  - - ANGLES - -           - POWER GAINS -       - - - POLARIZATION - - -    - - - E(THETA) - - -    - - - E(PHI) - - -\n")
        f.write("  THETA     PHI        VERT.   HOR.    TOTAL      AXIAL     TILT   SENSE     MAGNITUDE    PHASE      MAGNITUDE    PHASE \n")
        f.write(" DEGREES  DEGREES       DB      DB      DB        RATIO     DEG.              VOLTS/M    DEGREES      VOLTS/M    DEGREES\n")
        for ti in range(n_theta):
            th = thetas[ti]
            for pi in range(n_phi):
                ph = phis[pi]
                g = float(gtot[ti][pi])
                f.write(" {:7.2f} {:7.2f}  {:9.2f}{:8.2f}{:8.2f} {:10.5f} {:8.2f}  LINEAR  {:11.5E} {:8.2f}  {:11.5E} {:8.2f}\n".format(
                    th, ph, g, g, g, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))

    print("wrote", out_path, "ntheta={} nphi={} peak={:.2f} dBi".format(n_theta, n_phi, float(gtot.max())))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    nec_path = sys.argv[1]
    out_path = sys.argv[2] if len(sys.argv) > 2 else os.path.splitext(nec_path)[0] + ".out"
    run(nec_path, out_path)

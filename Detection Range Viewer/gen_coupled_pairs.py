"""
Bakes precomputed real-NEC mutual-coupling delta patterns for pairs of
vertically-stacked antennas, consumed by index.html's "Measured coupling"
combine mode (a third mode alongside "Union" and "Coherent (array-factor)").

Why a delta, and why free space:
Array-factor superposition (today's "Coherent" mode) phase-combines each
antenna's own independently-computed, unperturbed far-field pattern - it
can only reweight lobes that already exist in each isolated pattern, never
invent one. Real mutual coupling (the near-field presence of a second
antenna's actual wire geometry perturbing the fed antenna's own current
distribution) can grow large new lobes that don't exist in either antenna's
own isolated pattern at all - confirmed this session via PyNEC re-runs of
the user's own PLC1669_dual_opposing_*.nec files (isolated F/B ratio ~25 dB
collapses to ~6 dB at 1.0 m spacing, ~2.5-3.5 dB at 0.5 m).

Mutual coupling is a near-field, geometry-only effect between the two
physical structures; the app's live two-ray ground model
(computeLiveGrid/groundReflectionDb in index.html) is a separate far-field
effect that only depends on each antenna's own absolute height. So every
pair here is solved in free space (GE 0, no GN card) and reduced to
    deltaDb(elev,phi) = coupled_freeSpace_gainDb - isolated_freeSpace_gainDb
At runtime this delta is added on top of whatever live/baked grid the app
already shows for that instance's real height/terrain - the height slider
keeps working, and the delta isolates just the coupling contribution.

Geometry: each model's baseline single-antenna .nec source is parsed for its
GW (wire)/EX (feed)/FR (frequency) cards, then recentered so the boom's
geometric X-center sits at X=0 (the mast axis) - the same "mounted at their
midpoint" convention as the user's own PLC1669_dual_opposing_*.nec files
(shift so the boom center, not the driven element, sits above the mast).
Antenna B is that same recentered geometry rotated about the shared Z (mast)
axis by the pair's relative heading and shifted up by the spacing, with wire
tags offset by len(wires) to avoid collision (matching the user's own
dual_opposing files: tags 1-9 = antenna 1, 10-18 = antenna 2, tag 11 =
antenna 2's driven element). Only relative headings 0/45/90/180 deg are
solved - any other relative heading is the phi-mirror of 360-heading (a
left-right mirror of the whole two-antenna configuration, since both
baseline geometries are themselves mirror-symmetric about their own Z=const
boom plane), which index.html folds at lookup time instead of needing new
solves.

Both antennas are solved as the fed element in turn (role A = lower antenna
fed, upper parasite; role B = upper antenna fed, lower parasite) rather than
assumed interchangeable, producing deltaGridA/deltaGridB per (heading,
spacing) sample - one for whichever paired instance is lower, one for
whichever is upper.

Runs at a coarse 10 deg (theta) x 15 deg (phi) RP grid - the correction is a
smooth bump, not fine angular detail, so this keeps output size small;
index.html bilinearly interpolates it back up to each model's native
resolution. The 15 deg phi step (not 10) is deliberate: it evenly divides
every relative heading this script solves (45/90/180), which the rotation
correction below depends on. Uses this project's own PyNEC-based run_nec.py
to solve each deck; results are parsed with a small local parser (below)
rather than parse_nec.py's build_entry(), because build_entry() downsamples
phi to only even-integer degrees (phi_step=2), which would silently discard
every 15 deg-but-odd sample point (15, 45, 75, ...) - fine for the fine
native-resolution grids build_entry() was written for, but incompatible
with this script's coarser, heading-aligned grid.

IMPORTANT - role B rotation correction: antenna B's geometry is antenna A's
geometry rotated by the pair's relative heading about the shared mast (Z)
axis, so B's own natural (isolated) far-field pattern is A's isolated
pattern with phi shifted by +heading (rotating a radiator by theta rotates
its own pattern by theta, in the same lab-frame phi NEC reports). Diffing
role B's coupled output directly against iso_grid at the SAME lab phi (as
if B were unrotated) compares two different physical directions whenever
the isolated pattern has a sharp feature (e.g. a narrow null) - this was
caught via a smoke test that produced a spurious +30.7 dB "coupling" spike
at 1.5 lambda spacing (meant to be the near-zero-coupling anchor), traced to
a narrow -45.49 dBi null in the isolated pattern at phi=70 deg landing on a
grid point where role B's (differently-oriented) pattern has no comparable
null. The fix: deltaGridB is expressed in ANTENNA B's OWN local frame
(matching how instanceNativeGrid() already represents every instance - un-
rotated, with headingDeg applied later at resample time) via
    deltaGridB_local[Q] = roleB_labframe[(Q + heading) mod 360] - iso_grid[Q]
Since heading is always an exact multiple of the 15 deg phi step, this is
an exact index shift, never an interpolation. Role A needs no such
correction since antenna A is never rotated (its own local frame IS the lab
frame).

Usage: python gen_coupled_pairs.py [output.json]
If output.json is omitted, writes coupled_patterns.json next to this script.
Offline one-time bake: 2 models x 2 polarizations x (1 isolated + 4 headings
x 7 spacings x 2 roles) = 228 NEC solves total, each a tiny 10x24-point far
field (not the fine patterns.json grids), so this should finish in well
under a couple of minutes even though it's a lot of solves.
"""
import sys, os, json, math, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_nec
import parse_nec

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

BASE_FILES = {
    ("plc1669", "horizontal"): os.path.join(ROOT, "Laird PLC1669", "PLC1669_9el_166MHz.nec"),
    ("plc1669", "vertical"): os.path.join(ROOT, "Laird PLC1669", "PLC1669_9el_166MHz_vertical_freespace_grid.nec"),
    ("zda433", "horizontal"): os.path.join(ROOT, "ZDADJ433-12YG", "ZDADJ433_12YG_9el_433MHz_horizontal_freespace_grid.nec"),
    ("zda433", "vertical"): os.path.join(ROOT, "ZDADJ433-12YG", "ZDADJ433_12YG_9el_433MHz_vertical_freespace_grid.nec"),
}

HEADINGS_DEG = [0, 45, 90, 180]
SPACINGS_LAMBDA = [0.15, 0.25, 0.4, 0.6, 0.85, 1.15, 1.5]
RP_N_THETA, RP_DTHETA = 10, 10.0   # theta 0..90 step10 -> elevDeg 0..90 step10
RP_N_PHI, RP_DPHI = 25, 15.0        # phi 0..360 step15 (incl. dup 360) -> phiDeg 0..345 step15 after dedup
                                     # (15deg chosen so it evenly divides every
                                     # heading below - required for the exact
                                     # index-shift rotation correction on role B)


def load_base(path):
    """Extracts wire geometry + feed + frequency from a baseline single-
    antenna .nec, ignoring any GN/RP cards it has (this script builds its own
    free-space deck and RP grid). Recenters the boom so its geometric
    X-center sits at X=0 (the shared mast axis), matching the "mounted at
    their midpoint" convention in the user's own dual-antenna example
    files."""
    comments, cards = run_nec.parse_nec_file(path)
    wires = []
    fed_tag = fed_seg = freq_mhz = None
    for tag, toks in cards:
        nums = [float(t) for t in toks]
        if tag == "GW":
            wires.append({
                "segs": int(nums[1]),
                "x1": nums[2], "y1": nums[3], "z1": nums[4],
                "x2": nums[5], "y2": nums[6], "z2": nums[7],
                "rad": nums[8],
            })
        elif tag == "EX":
            fed_tag, fed_seg = int(nums[1]), int(nums[2])
        elif tag == "FR":
            freq_mhz = nums[4]
    assert wires and fed_tag is not None and freq_mhz, "couldn't parse geometry/feed/freq from " + path
    xs = [w["x1"] for w in wires] + [w["x2"] for w in wires]
    xshift = -(min(xs) + max(xs)) / 2.0
    for w in wires:
        w["x1"] += xshift
        w["x2"] += xshift
    return wires, fed_tag, fed_seg, freq_mhz


def transform_wires(wires, dz=0.0, heading_deg=0.0):
    """Rotates each wire about the shared Z (mast) axis by heading_deg, then
    shifts it up by dz - used to place antenna B relative to antenna A,
    which is left untouched at the recentered baseline geometry."""
    th = math.radians(heading_deg)
    cos_t, sin_t = math.cos(th), math.sin(th)

    def rot(x, y):
        return x * cos_t - y * sin_t, x * sin_t + y * cos_t

    out = []
    for w in wires:
        x1, y1 = rot(w["x1"], w["y1"])
        x2, y2 = rot(w["x2"], w["y2"])
        out.append({
            "segs": w["segs"],
            "x1": x1, "y1": y1, "z1": w["z1"] + dz,
            "x2": x2, "y2": y2, "z2": w["z2"] + dz,
            "rad": w["rad"],
        })
    return out


def write_deck(path, wires_a, wires_b, fed_tag, fed_seg, freq_mhz, comment):
    lines = ["CM " + comment, "CE"]
    tag = 1
    for w in wires_a + wires_b:
        lines.append("GW {} {} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.7f}".format(
            tag, w["segs"], w["x1"], w["y1"], w["z1"], w["x2"], w["y2"], w["z2"], w["rad"]))
        tag += 1
    lines.append("GE 0")
    lines.append("EX 0 {} {} 0 1 0".format(fed_tag, fed_seg))
    lines.append("FR 0 1 0 0 {:.4f} 0".format(freq_mhz))
    lines.append("RP 0 {} {} 1000 0.0 0.0 {:.1f} {:.1f} 0.0 0.0".format(RP_N_THETA, RP_N_PHI, RP_DTHETA, RP_DPHI))
    lines.append("EN")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def parse_coarse_out(path):
    """Parses a .out from write_deck's fixed RP grid (theta 0..90 step
    RP_DTHETA, phi 0..360 step RP_DPHI) into {elevDeg, phiDeg, gainDb}.

    Deliberately NOT parse_nec.build_entry(): that function's phi_step=2
    downsampling assumes a fine native grid and would silently drop every
    odd-degree phi sample (15, 45, 75, ...) that this script's 15deg-step
    grid depends on for the role-B rotation correction. Also skips
    build_entry()'s zenith auto-detection since this script's own RP card
    is known in advance to sweep theta=0 (zenith) -> 90 (horizon)."""
    with open(path, "r", errors="ignore") as f:
        lines = f.readlines()
    rows = parse_nec.parse_pattern_rows(lines)
    assert rows, "no RADIATION PATTERNS section in " + path

    thetas = sorted(set(round(t, 2) for t, p, g in rows))
    phis_all = sorted(set(round(p, 2) for t, p, g in rows))
    phis = [p for p in phis_all if p < 360.0] if (phis_all[0] == 0.0 and phis_all[-1] == 360.0) else phis_all

    grid = {}
    for t, p, g in rows:
        grid.setdefault((round(t, 2), round(p % 360.0, 2)), g)

    def gain_at(t, p):
        g = grid.get((round(t, 2), round(p % 360.0, 2)))
        return -60.0 if g is None else g

    def sanitized(t, p):
        g = gain_at(t, p)
        if g > -900.0:
            return g
        # NEC-2 emits a -999.99 sentinel at a handful of exact-null grid
        # points rather than a real value - fall back to the average of the
        # same theta's immediate phi neighbors (each phi step apart).
        step = phis[1] - phis[0] if len(phis) > 1 else RP_DPHI
        neighbors = [gain_at(t, p - step), gain_at(t, p + step)]
        neighbors = [v for v in neighbors if v > -900.0]
        return sum(neighbors) / len(neighbors) if neighbors else g

    elevs = sorted(set(round(90.0 - t, 2) for t in thetas))
    matrix = [[round(sanitized(90.0 - e, p), 2) for p in phis] for e in elevs]
    return {"elevDeg": elevs, "phiDeg": phis, "gainDb": matrix}


def solve(deck_writer_args, tmpdir, tag):
    nec_path = os.path.join(tmpdir, tag + ".nec")
    out_path = os.path.join(tmpdir, tag + ".out")
    write_deck(nec_path, *deck_writer_args)
    run_nec.run(nec_path, out_path)
    return parse_coarse_out(out_path)


def main():
    out_path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "coupled_patterns.json")
    tmpdir = tempfile.mkdtemp(prefix="coupled_nec_")
    print("scratch dir:", tmpdir)

    result = {}
    for model, pol in BASE_FILES:
        wires, fed_tag, fed_seg, freq_mhz = load_base(BASE_FILES[(model, pol)])
        lambda_m = 299.792458 / freq_mhz
        print("== {} / {} ({:.3f} MHz, lambda={:.3f} m, {} wires, fed tag {}) ==".format(
            model, pol, freq_mhz, lambda_m, len(wires), fed_tag))

        iso_entry = solve((wires, [], fed_tag, fed_seg, freq_mhz,
                                  "isolated {} {} baseline (coarse grid, coupling-delta reference)".format(model, pol)),
                           tmpdir, "{}_{}_isolated".format(model, pol))
        coarse_elev, coarse_phi = iso_entry["elevDeg"], iso_entry["phiDeg"]
        iso_grid = iso_entry["gainDb"]
        nE, nP = len(coarse_elev), len(coarse_phi)

        pol_out = {"coarseElevDeg": coarse_elev, "coarsePhiDeg": coarse_phi, "headings": {}}
        for heading in HEADINGS_DEG:
            samples = []
            for spacing_lambda in SPACINGS_LAMBDA:
                spacing_m = spacing_lambda * lambda_m
                wires_b = transform_wires(wires, dz=spacing_m, heading_deg=heading)
                tag_suffix = "{}_{}_{}_{}".format(model, pol, heading, spacing_lambda)

                # At heading=0 (directly stacked, no azimuth offset), a small
                # enough spacing can put antenna B's elements at the exact
                # same X/Y as antenna A's, overlapping in Z once the spacing
                # is shorter than an element's own Z-extent (only reachable
                # for "vertical" polarization, whose elements run along Z) -
                # a genuine physical impossibility for this simple co-axial
                # stacking geometry, not a bug. PyNEC raises on the resulting
                # intersecting wires; skip that one sample rather than fail
                # the whole bake. index.html's nearest-neighbor snap will
                # fall back to the closest spacing that DID solve.
                try:
                    fed_tag_b = fed_tag + len(wires)
                    entry_a = solve((wires, wires_b, fed_tag, fed_seg, freq_mhz,
                                            "coupled pair role=A(lower fed) model={} pol={} heading={} spacingLambda={}".format(
                                                model, pol, heading, spacing_lambda)),
                                     tmpdir, "roleA_" + tag_suffix)
                    entry_b = solve((wires, wires_b, fed_tag_b, fed_seg, freq_mhz,
                                            "coupled pair role=B(upper fed) model={} pol={} heading={} spacingLambda={}".format(
                                                model, pol, heading, spacing_lambda)),
                                     tmpdir, "roleB_" + tag_suffix)
                except RuntimeError as e:
                    print("  heading={:>3} spacing={:.2f}lambda ({:.2f}m)  SKIPPED (geometry infeasible: {})".format(
                        heading, spacing_lambda, spacing_m, e))
                    continue

                assert entry_a["elevDeg"] == coarse_elev and entry_a["phiDeg"] == coarse_phi
                assert entry_b["elevDeg"] == coarse_elev and entry_b["phiDeg"] == coarse_phi

                # Antenna B's geometry is antenna A's rotated by `heading` deg
                # about the shared mast axis, so B's own isolated pattern would
                # be iso_grid shifted by +heading in lab-frame phi. Comparing
                # roleB's coupled output straight against iso_grid at the same
                # lab phi (see module docstring) mismatches physical directions
                # whenever iso_grid has a sharp feature. Un-rotate roleB's
                # output back into its own local frame before differencing:
                # deltaGridB_local[Q] = roleB_lab[(Q+heading)] - iso_grid[Q].
                # heading is always an exact multiple of the phi step, so this
                # is an exact index shift, never an interpolation.
                phi_step = coarse_phi[1] - coarse_phi[0]
                shift = int(round(heading / phi_step)) % nP
                assert abs(shift * phi_step - heading) < 1e-6, "heading not a multiple of the phi grid step"

                delta_a = [[round(entry_a["gainDb"][ti][pi] - iso_grid[ti][pi], 2) for pi in range(nP)] for ti in range(nE)]
                delta_b = [[round(entry_b["gainDb"][ti][(pi + shift) % nP] - iso_grid[ti][pi], 2) for pi in range(nP)] for ti in range(nE)]
                peak_a = max(max(row) for row in delta_a)
                peak_b = max(max(row) for row in delta_b)
                print("  heading={:>3} spacing={:.2f}lambda ({:.2f}m)  max deltaA={:+.2f}dB max deltaB={:+.2f}dB".format(
                    heading, spacing_lambda, spacing_m, peak_a, peak_b))
                samples.append({"spacingLambda": spacing_lambda, "deltaGridA": delta_a, "deltaGridB": delta_b})
            pol_out["headings"][str(heading)] = samples
        result.setdefault(model, {})[pol] = pol_out

    with open(out_path, "w") as f:
        json.dump(result, f)
    print("wrote", out_path, "(" + str(os.path.getsize(out_path)) + " bytes)")


if __name__ == "__main__":
    main()

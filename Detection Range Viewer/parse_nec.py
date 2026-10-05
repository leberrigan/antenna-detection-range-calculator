"""
Convert one or more folders of NEC-2 .out radiation-pattern files into
patterns.json for the 3D detection-range viewer (index.html) in this same
folder. Each folder is a distinct antenna model (e.g. a 166 MHz 9-el Yagi
vs. a synthesized 434 MHz 7-el Yagi) and becomes its own top-level entry,
selectable in the viewer's model dropdown.

Usage:
    python parse_nec.py "<folder 1>" ["<folder 2>" ...] [output.json]

If output.json is omitted, patterns.json is written next to this script.

What it does with each folder:
  - Assigns it a model key + display label from MODEL_LABELS below (matched
    by folder name), falling back to a slug/title of the folder name for an
    unlisted folder.

What it does with each *.out file within a folder:
  - A file with no GN (ground plane) card is treated as the free-space
    baseline pattern (key "freeSpace") that the viewer's live height/terrain
    model builds on, rather than a terrain - but only if it's a full
    theta/phi grid (same RP card shape as the terrain runs); a partial
    sweep (e.g. a two-cut E/H-plane sanity check) is still skipped.
  - Reads the ground dielectric constant (eps) and conductivity (sigma) from
    the echoed "DATA CARD ... GN ..." line.
  - Reads a friendly ground label from a "GROUND TYPE: <label> (...)" comment
    if the .nec file has one (CM cards echo into the .out); otherwise falls
    back to eps/sigma, or the filename.
  - Auto-detects which end of the theta sweep is the zenith (straight up):
    that's the end where gain barely changes with phi, since straight up is
    a single physical direction. The other end is the horizon. This makes
    the elevation convention self-correcting regardless of how a given NEC
    run happened to sign/offset its THETA values.
  - Re-expresses every row by elevation in degrees, 0 = horizon, 90 = zenith,
    sorted ascending, so the viewer needs no per-dataset convention logic.
  - Detects polarization (horizontal/vertical) from the filename or a "...
    polarization" comment, and groups entries in the output as
    {polarization: {terrainKey: entry, ...}, ...} so the viewer can toggle
    between polarizations while keeping the same terrain selected. Terrain
    keys are assigned by matching each file's ground eps/sigma against a
    table of known terrains, so the same physical ground gets the same key
    regardless of which polarization's file it came from.

Re-run this any time you add/change .nec/.out files for an antenna, then
reopen index.html (or republish it) to see the updated patterns.
"""
import re, json, os, sys, glob

FLOAT_RE = re.compile(r'^-?\d+\.?\d*(?:[Ee][+-]?\d+)?$')

# known terrains by (eps, sigma), so files from different polarization runs
# that describe the same ground collapse to the same key in the output
CANONICAL_TERRAINS = [
    ("average",    13, 0.005),
    ("desert",     3,  0.00015),
    ("forest",     13, 0.006),
    ("freshwater", 80, 0.01),
    ("grassland",  14, 0.01),
    ("seawater",   81, 5.0),
]

# known antenna-model folders, keyed by folder basename (case-insensitive) ->
# (model key used in patterns.json, display label shown in the viewer's
# dropdown). An unlisted folder falls back to a slug/title of its own name.
MODEL_LABELS = {
    "laird plc1669":          ("plc1669",  "PLC1669 · 9-EL YAGI"),
    "zdadj433-12yg":          ("zda433",   "ZDA ZDADJ433-12YG · 9-EL YAGI"),
    "wa5vjb cheapyagi":       ("wa5vjb_yagi", "WA5VJB CHEAP YAGI · 11-EL YAGI (2.45GHz)"),
    "hg2412p corner reflector": ("hg2412p", "HG2412P · CORNER REFLECTOR (2.45GHz)"),
    "zdaqj166":                ("zdaqj166", "ZDAQJ166 · OMNI DIPOLE (166MHz)"),
    "dx 9el yagi":             ("dx_9el_yagi", "DX SHOP 9-EL YAGI (166MHz)"),
}


def model_info(folder):
    name = os.path.basename(os.path.normpath(folder))
    hit = MODEL_LABELS.get(name.lower())
    if hit:
        return hit
    slug = re.sub(r'[^a-z0-9]+', '_', name.lower()).strip('_') or name
    return slug, name.title()


def canonical_terrain_key(ground):
    if not ground:
        return None
    for slug, eps, sigma in CANONICAL_TERRAINS:
        sigma_tol = max(0.0004, sigma * 0.05)
        if abs(ground["eps"] - eps) < 0.5 and abs(ground["sigma"] - sigma) < sigma_tol:
            return slug
    return None


def detect_polarization(path, lines):
    stem = os.path.splitext(os.path.basename(path))[0]
    tokens = [t.lower() for t in re.split(r'[^A-Za-z0-9]+', stem) if t]
    text = " ".join(lines)
    if "cp" in tokens or re.search(r'circular\s*polarization', text, re.I):
        return "circular"
    if "vertical" in tokens or re.search(r'vertical\s*polarization', text, re.I):
        return "vertical"
    return "horizontal"


def try_parse_row(line):
    toks = line.split()
    if len(toks) < 11:
        return None
    head = toks[:7]
    if not all(FLOAT_RE.match(t) for t in head):
        return None
    theta, phi, vert, hor, total, axial, tilt = (float(t) for t in head)
    return theta, phi, total


def parse_pattern_rows(lines):
    start = None
    for i, l in enumerate(lines):
        if 'RADIATION PATTERNS' in l:
            start = i
            break
    if start is None:
        return None
    rows = []
    started = False
    for l in lines[start + 1:]:
        parsed = try_parse_row(l)
        if parsed is not None:
            started = True
            rows.append(parsed)
        elif started:
            if l.strip() == '':
                continue
            break
    return rows


def parse_ground(lines):
    for l in lines:
        if 'DATA CARD' in l and re.search(r'\bGN\b', l):
            toks = l.split()
            try:
                gi = toks.index('GN')
                eps = float(toks[gi + 5])
                sigma = float(toks[gi + 6])
                return {"eps": eps, "sigma": sigma}
            except (ValueError, IndexError):
                return None
    return None


def parse_ground_label(lines):
    for l in lines:
        m = re.search(r'GROUND TYPE:\s*(.+)', l, re.I)
        if m:
            # drop the parenthetical detail - some CM comment lines get cut off
            # mid-parenthetical by NEC's line-length limit, so don't require it
            # to actually close
            label = re.sub(r'\s*\(.*$', '', m.group(1)).strip()
            if label == label.upper():
                label = label.title()
            return label
    return None


def parse_freq_mhz(lines):
    for l in lines:
        m = re.search(r'FREQUENCY\s*=\s*([\d.Ee+-]+)\s*MHZ', l, re.I)
        if m:
            return float(m.group(1))
    return None


def default_label(path, ground):
    if ground and abs(ground["eps"] - 13) < 0.5 and abs(ground["sigma"] - 0.005) < 0.001:
        return "Average Ground"
    stem = os.path.splitext(os.path.basename(path))[0]
    words = re.sub(r'[_\-]+', ' ', stem).strip()
    return words.title() if words else stem


def build_entry(path):
    with open(path, 'r', errors='ignore') as f:
        lines = f.readlines()

    # ground is None for a free-space run (no GN card) - that's not skipped
    # here anymore, since the viewer's analytic height/terrain model needs a
    # full-grid free-space baseline pattern as its starting point. A GN-less
    # file that *isn't* a full grid (e.g. a two-cut E/H-plane sanity check)
    # is still rejected below once the grid shape is known.
    ground = parse_ground(lines)

    rows = parse_pattern_rows(lines)
    if not rows:
        return None, "no RADIATION PATTERNS section found"

    # NEC-2's theta is only ever physically in [0,180] (angle from +Z), but a
    # custom RP sweep can use a negative theta0 (e.g. -90..0 or -90..+90) to
    # cover the same hemisphere from "below" instead of extending past it.
    # A row recorded at negative theta t is the exact same physical direction
    # as (|t|, phi+180): sin(t)=-sin(|t|) while cos(t)=cos(|t|), which flips
    # the sign of both x and y (i.e. rotates phi by 180deg) without changing
    # z. Left uncorrected, a file whose sweep never crosses back to positive
    # theta (so the dedup below never discards its negative half) keeps every
    # row's raw phi as-is, mirroring its whole azimuth pattern by 180deg.
    rows = [(t, ((p + 180.0) % 360.0) if t < 0 else p, g) for (t, p, g) in rows]

    thetas = sorted(set(r[0] for r in rows))
    phis = sorted(set(r[1] for r in rows))
    phis_used = [p for p in phis if p < 360.0] if (phis[0] == 0.0 and phis[-1] == 360.0) else phis
    if len(thetas) < 5 or len(phis_used) < 5:
        return None, "not a full theta/phi pattern grid (e.g. a single E/H-plane cut)"
    phi_step = 2
    phis_ds = [p for p in phis_used if int(round(p)) % phi_step == 0] or phis_used

    grid = {(round(t, 2), round(p, 2)): g for (t, p, g) in rows}

    def gain_at(t, p):
        g = grid.get((round(t, 2), round(p, 2)))
        if g is None:
            g = grid.get((round(t, 2), round(p % 360.0, 2)))
        return g if g is not None else -60.0

    # auto-detect the zenith of the theta sweep: straight up is one physical
    # direction, so gain there is ~constant across phi, unlike the horizon
    # end which varies strongly with phi (compass direction). Some sweeps
    # run past zenith on both sides (e.g. theta -90..+90 instead of -90..0):
    # NEC-2's negative-theta convention makes the far side an exact mirror
    # of the near side (gain(t,p) == gain(-t,p+180)), i.e. the same upper
    # hemisphere counted twice rather than new below-ground data, so only
    # the half-sweep from the zenith to its farther end is kept.
    def phi_spread(t):
        vals = [gain_at(t, p) for p in phis_ds]
        # NEC-2 emits a -999.99 sentinel for a handful of exact-null grid
        # points (e.g. theta=0/phi=0 on some geometries) rather than a real
        # deep-fade gain value. A single one of these dwarfs any genuine
        # phi-to-phi variation and fools the zenith pick below into treating
        # a real horizon/off-axis theta as if it were the phi-invariant pole,
        # so it's excluded here the same way a missing grid point already is.
        real_vals = [v for v in vals if v > -900.0]
        vals = real_vals or vals
        return max(vals) - min(vals)

    # NEC-2's theta is defined as angle from +Z by fixed convention (not a
    # per-file choice), and every model in this project is built with Z=up,
    # so theta=0.0 is always the true zenith when it's actually in the swept
    # grid - no need to infer it from phi variation, which can be genuinely
    # noisy/large near a real pattern null on-axis (e.g. a Z-oriented dipole
    # element has a deep null straight up, its own antenna axis) rather than
    # small the way it is near a lobe. Only fall back to the phi-spread guess
    # when 0.0 isn't itself one of the swept theta samples.
    zenith_theta = 0.0 if 0.0 in thetas else min(thetas, key=phi_spread)
    lo, hi = thetas[0], thetas[-1]
    if zenith_theta not in (lo, hi):
        thetas = [t for t in thetas if t >= zenith_theta] if (hi - zenith_theta) >= (zenith_theta - lo) \
            else [t for t in thetas if t <= zenith_theta]
    # the branch above only dedupes a sweep that runs past zenith on *both*
    # sides (zenith strictly interior); a sweep that runs past zenith on only
    # one side but by more than 90 deg (e.g. -180..0, with zenith already
    # sitting at the far endpoint) isn't caught by it and would otherwise
    # leave elevDeg going negative below - a real one-hemisphere-only theta
    # range never needs more than 90 deg of "distance" from zenith, so this
    # clamp is a safe no-op for every already-correct dataset.
    thetas = [t for t in thetas if abs(t - zenith_theta) <= 90]

    elev_rows = sorted(
        ((round(90 - abs(t - zenith_theta), 3), t) for t in thetas),
        key=lambda x: x[0]
    )

    def sanitized_gain_at(t, p):
        # Same -999.99 NEC-2 sentinel as phi_spread() above, but here it would
        # otherwise land directly in the exported pattern as a bogus deep
        # gouge (skewing peak/min stats and the Option A free-space
        # baseline). Patch it with the average of the other phi cuts at the
        # same theta, which is a good estimate since real gain varies
        # smoothly with phi and is near-constant right at the pole anyway.
        g = gain_at(t, p)
        if g <= -900.0:
            others = [gain_at(t, pp) for pp in phis_ds if pp != p]
            others = [v for v in others if v > -900.0]
            if others:
                return sum(others) / len(others)
        return g

    elevs = [er[0] for er in elev_rows]
    matrix = [[round(sanitized_gain_at(t, p), 2) for p in phis_ds] for (_, t) in elev_rows]

    label = parse_ground_label(lines) or default_label(path, ground)
    return {
        "label": label,
        "ground": ground,
        "freqMhz": parse_freq_mhz(lines),
        "polarization": detect_polarization(path, lines),
        "elevDeg": elevs,     # ascending: 0 = horizon .. 90 = zenith
        "phiDeg": phis_ds,
        "gainDb": matrix,
    }, None


def build_manifest_for_folder(folder):
    """Returns {polarization: {terrainKey: entry}} for one antenna model's folder."""
    pol_manifest = {}
    for path in sorted(glob.glob(os.path.join(folder, "*.out"))):
        stem = os.path.splitext(os.path.basename(path))[0]
        entry, skip_reason = build_entry(path)
        if entry is None:
            print("skip {}: {}".format(os.path.basename(path), skip_reason))
            continue
        pol = entry.pop("polarization")
        if entry["ground"] is None:
            key = "freeSpace"
        else:
            key = canonical_terrain_key(entry["ground"]) or re.sub(r'[^a-z0-9]+', '_', stem.lower()).strip('_') or stem
        pol_manifest.setdefault(pol, {})[key] = entry
        print("{:10s} {:24s} -> {:28s}  {} theta x {} phi  peak {:.1f} dBi".format(
            pol, key, entry["label"], len(entry["elevDeg"]), len(entry["phiDeg"]),
            max(max(row) for row in entry["gainDb"])
        ))
    return pol_manifest


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    args = sys.argv[1:]
    if args[-1].lower().endswith(".json"):
        out_path = args[-1]
        folders = args[:-1]
    else:
        out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "patterns.json")
        folders = args

    manifest = {}
    for folder in folders:
        model_key, model_label = model_info(folder)
        print("== {} ({}) ==".format(model_label, model_key))
        pol_manifest = build_manifest_for_folder(folder)
        if not pol_manifest:
            print("  no usable ground-plane radiation patterns found, skipping model")
            continue
        manifest[model_key] = {"label": model_label, "polarizations": pol_manifest}

    if not manifest:
        print("No usable ground-plane radiation patterns found in", folders)
        sys.exit(1)

    with open(out_path, "w") as f:
        json.dump(manifest, f)
    print("wrote", out_path, "(" + str(os.path.getsize(out_path)) + " bytes)")


if __name__ == "__main__":
    main()

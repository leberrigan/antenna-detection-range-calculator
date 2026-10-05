"""
Local HTTP server backing index.html's "Live coupling" combine mode: solves
real mutual coupling for the EXACT current antenna configuration (arbitrary
heights/headings, any N from 2 up to MAX_ANTENNAS, even mixing models that
share a frequency) on demand via PyNEC, instead of looking up the
precomputed 2-antenna/discretized-heading library in coupled_patterns.json
(gen_coupled_pairs.py / "Measured coupling" mode - untouched by this file).

Also serves index.html and its sibling static files (patterns.json,
coupled_patterns.json, ...) from the same origin, so this single process
replaces `python -m http.server` as the way to open the app - no CORS setup
needed. Opening index.html via file:// or plain `python -m http.server`
still works for everything except "Live coupling", which the page disables
with an explanatory note when this server (or its /compute-coupling
endpoint) isn't reachable.

Usage: python coupling_server.py [port]   (default 8765)
Run from this folder (Detection Range Viewer/); it chdir()s here anyway so
it also works if invoked from elsewhere.

Physics / math, in brief (see gen_coupled_pairs.py's module docstring for
the full background on why a free-space delta, and the role-B rotation
problem this generalizes):
  - Every participating instance's baseline wire geometry (one raw .nec per
    (model, polarization), BASELINE_FILES below) is loaded once via
    gen_coupled_pairs.load_base(), which recenters its boom to X=0 (mast
    axis) - untouched, reused as-is.
  - Each instance is placed in one shared lab frame via
    gen_coupled_pairs.transform_wires(wires, dz, heading_deg), using the
    instance's REAL absolute heightM/headingDeg directly (not relative to a
    reference antenna, unlike the offline pairwise script, which only ever
    needed relative spacing/heading between exactly 2 same-model antennas).
    dz is heightM minus that baseline file's own nominal Z (mean wire Z),
    so different models' baselines - each possibly baked at a slightly
    different nominal height - align onto one real-world height axis.
  - One joint NEC deck contains every participant's wires. NEC can only feed
    one segment per solve, so a full N-antenna joint solve needs N solves
    (one per instance as the fed/"role" element, every other instance
    present as passive parasitic geometry) - generalizing
    gen_coupled_pairs.py's role-A/role-B pattern from 2 to N.
  - Rotation trick: request each role's RP card with phi0 = that instance's
    own heading_deg (mod 360), instead of always phi0=0 like the offline
    script. NEC samples phi = phi0 + k*dphi, i.e. physical/lab angle
    heading_i + k*dphi at index k. Since instance i's geometry is the
    baseline rotated by heading_i, its own local-frame angle Q maps to lab
    angle Q+heading_i - so index k (lab angle heading_i+k*dphi) is exactly
    local-frame angle k*dphi, for ANY real heading_i, with no post-hoc
    index-shift or interpolation (unlike the offline script, which could
    only afford an exact index shift because heading was restricted to
    multiples of the 15 deg phi step). This is what lets heights/headings
    stay fully continuous in this mode.
  - Each instance's delta = its role-solve output minus the ISOLATED
    baseline for its own (model, polarization) - solved once (phi0=0, dz=0,
    no other geometry present) and cached in-process per (model, pol),
    reused across requests and across every instance that shares a model.
"""
import sys, os, json, math, itertools, traceback
import http.server

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen_coupled_pairs as gcp
from PyNEC import nec_context

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

BASELINE_FILES = {
    ("plc1669", "horizontal"): os.path.join(ROOT, "Laird PLC1669", "PLC1669_9el_166MHz.nec"),
    ("plc1669", "vertical"): os.path.join(ROOT, "Laird PLC1669", "PLC1669_9el_166MHz_vertical_freespace_grid.nec"),
    ("zda433", "horizontal"): os.path.join(ROOT, "ZDADJ433-12YG", "ZDADJ433_12YG_9el_433MHz_horizontal_freespace_grid.nec"),
    ("zda433", "vertical"): os.path.join(ROOT, "ZDADJ433-12YG", "ZDADJ433_12YG_9el_433MHz_vertical_freespace_grid.nec"),
    ("wa5vjb_yagi", "horizontal"): os.path.join(ROOT, "WA5VJB CheapYagi", "WA5VJB_CheapYagi_11el_2450MHz_horizontal_freespace_grid.nec"),
    ("hg2412p", "vertical"): os.path.join(ROOT, "HG2412P Corner Reflector", "HG2412P_derived_corner_reflector_2450MHz_vertical_freespace_grid.nec"),
    ("zdaqj166", "vertical"): os.path.join(ROOT, "ZDAQJ166", "ZDAQJ166_166MHz_vertical_freespace_grid.nec"),
    ("dx_9el_yagi", "horizontal"): os.path.join(ROOT, "DX 9EL YAGI", "DX_9EL_YAGI_166MHz_freespace_grid.nec"),
}

RP_N_THETA, RP_DTHETA = 10, 10.0   # matches gen_coupled_pairs.py's coarse grid,
RP_N_PHI, RP_DPHI = 24, 15.0       # so index.html's bilinearCoupledDelta() can be
                                    # reused unmodified (no dup-360 point needed
                                    # here since nothing parses a .out file).
FREQ_TOL_MHZ = 0.05

_baseline_cache = {}   # (model, pol) -> (wires, fed_tag, fed_seg, freq_mhz, nominal_z, butt_x)
_iso_cache = {}         # (model, pol) -> {"elevDeg","phiDeg","gainDb"}


class SolveError(Exception):
    def __init__(self, code, message, detail=None):
        self.code, self.message, self.detail = code, message, detail


def get_baseline(model, pol):
    key = (model, pol)
    if key in _baseline_cache:
        return _baseline_cache[key]
    path = BASELINE_FILES.get(key)
    if not path or not os.path.isfile(path):
        raise SolveError("missing_baseline", "No live-solve wire geometry available for {} / {}".format(model, pol),
                          {"model": model, "polarization": pol})
    wires, fed_tag, fed_seg, freq_mhz, butt_x = gcp.load_base(path)
    zs = [w["z1"] for w in wires] + [w["z2"] for w in wires]
    nominal_z = sum(zs) / len(zs)
    entry = (wires, fed_tag, fed_seg, freq_mhz, nominal_z, butt_x)
    _baseline_cache[key] = entry
    return entry


def shift_wires_x(wires, dx):
    """Translates every wire by dx along X only (no rotation/Z change) -
    used to re-reference a "butt mounted" instance's boom so its butt end
    (rather than its geometric midpoint) sits on the shared mast axis,
    before transform_wires() rotates it into place for a given heading."""
    return [dict(w, x1=w["x1"] + dx, x2=w["x2"] + dx) for w in wires]


def solve_pattern(wires, fed_tag, fed_seg, freq_mhz, phi0):
    """In-memory equivalent of run_nec.run()'s PyNEC body (free-space deck,
    single EX/FR/RP card set) - skips the .nec/.out disk round-trip
    entirely. Returns (elevDeg[], phiDeg[], gainDb[][]), elevDeg ascending
    (0=horizon .. 90=zenith, matching run_nec.py/patterns.json convention).
    """
    try:
        c = nec_context()
        geo = c.get_geometry()
        for i, w in enumerate(wires, start=1):
            geo.wire(i, w["segs"], w["x1"], w["y1"], w["z1"], w["x2"], w["y2"], w["z2"], w["rad"], 1.0, 1.0)
        c.geometry_complete(0)   # free space, no ground card
        c.ex_card(0, fed_tag, fed_seg, 0, 1, 0, 0, 0, 0, 0)
        c.fr_card(0, 1, freq_mhz, 0)
        c.rp_card(0, RP_N_THETA, RP_N_PHI, 0, 0, 0, 0, 0.0, phi0, RP_DTHETA, RP_DPHI, 1000, 0.0)
        rp = c.get_radiation_pattern(0)
        thetas = rp.get_theta_angles()
        gtot = rp.get_gain_tot()
    except Exception as e:
        raise SolveError("infeasible_geometry", "NEC solve failed (likely intersecting wire geometry): " + str(e))

    import numpy as np
    gtot = np.array(gtot)
    if gtot.ndim == 1:
        # same phi-major/theta-minor reshape gotcha as run_nec.py - see its
        # comment for how this was verified against a live PyNEC run.
        gtot = gtot.reshape(RP_N_PHI, RP_N_THETA).T

    elev_deg = [round(90.0 - t, 2) for t in thetas]
    phi_deg = [round(k * RP_DPHI, 2) for k in range(RP_N_PHI)]
    raw = [[round(float(gtot[ti][pi]), 2) for pi in range(RP_N_PHI)] for ti in range(RP_N_THETA)]

    # NEC-2 emits a -999.99 sentinel at a handful of exact-null grid points
    # rather than a real value (same behavior gen_coupled_pairs.py's
    # sanitized() works around) - fall back to averaging the same theta's
    # immediate phi neighbors (grid wraps at 360, hence the modulo).
    def sanitized(ti, pi):
        g = raw[ti][pi]
        if g > -900.0:
            return g
        neighbors = [raw[ti][(pi - 1) % RP_N_PHI], raw[ti][(pi + 1) % RP_N_PHI]]
        neighbors = [v for v in neighbors if v > -900.0]
        if neighbors:
            return round(sum(neighbors) / len(neighbors), 2)
        # Whole phi ring degenerate - only happens at the theta=0 zenith
        # pole, a coordinate singularity where every phi maps to the same
        # physical direction (confirmed against patterns.json, which shows
        # the same all-sentinel row for e.g. hg2412p/vertical at elev=90).
        # Falling back to leaving -999.99 in place would only cancel
        # correctly if BOTH sides of a later delta subtraction hit the pole
        # identically - true for an isolated solve compared against itself,
        # but not once mutual coupling perturbs the joint solve enough to
        # break the isolated antenna's exact symmetry there. Borrow the
        # nearest elevation ring's average instead so the pole always has a
        # physically reasonable (if approximate) value.
        for other_ti in sorted(range(RP_N_THETA), key=lambda t: abs(t - ti)):
            if other_ti == ti:
                continue
            ring = [v for v in raw[other_ti] if v > -900.0]
            if ring:
                return round(sum(ring) / len(ring), 2)
        return g

    grid = [[sanitized(ti, pi) for pi in range(RP_N_PHI)] for ti in range(RP_N_THETA)]

    # thetas (and therefore elev_deg/grid rows, both indexed by ti) come out
    # of NEC in theta-ascending order, i.e. elev_deg descending (90=zenith
    # first). Reverse both so elev_deg is ascending (0=horizon..90=zenith),
    # matching this function's own documented contract and the
    # gen_coupled_pairs.py/coupled_patterns.json convention that index.html's
    # bilinearCoupledDelta() assumes - it clamps targetElev via
    # Math.max(coarseElev[0], Math.min(coarseElev[last], targetElev)), which
    # silently collapses every target to a constant elev when fed a
    # descending array instead of interpolating properly.
    elev_deg = list(reversed(elev_deg))
    grid = list(reversed(grid))
    return elev_deg, phi_deg, grid


def get_isolated(model, pol):
    key = (model, pol)
    if key in _iso_cache:
        return _iso_cache[key]
    wires, fed_tag, fed_seg, freq_mhz, _nominal_z, _butt_x = get_baseline(model, pol)
    elev_deg, phi_deg, grid = solve_pattern(wires, fed_tag, fed_seg, freq_mhz, 0.0)
    entry = {"elevDeg": elev_deg, "phiDeg": phi_deg, "gainDb": grid}
    _iso_cache[key] = entry
    return entry


def compute_coupling(antennas):
    """antennas: list of {"id","model","polarization","heightM","headingDeg"}.
    Returns the success response dict, or raises SolveError."""
    if len(antennas) < 2:
        raise SolveError("too_few_antennas", "Need at least 2 antennas to compute coupling")

    missing = [a["id"] for a in antennas if (a["model"], a["polarization"]) not in BASELINE_FILES]
    if missing:
        raise SolveError("missing_baseline", "No live-solve wire geometry for these instances", {"ids": missing})

    baselines = {a["id"]: get_baseline(a["model"], a["polarization"]) for a in antennas}
    freqs = {a["id"]: baselines[a["id"]][3] for a in antennas}
    ref_freq = next(iter(freqs.values()))
    mismatched = {aid: f for aid, f in freqs.items() if abs(f - ref_freq) > FREQ_TOL_MHZ}
    if mismatched:
        raise SolveError("frequency_mismatch", "All antennas in one live solve must share a frequency",
                          {"freqMhz": freqs})

    # place every participant in one shared lab frame, contiguous global wire
    # tags, tracking each instance's own global fed tag
    participants = []
    offset = 0
    for a in antennas:
        wires, fed_tag, fed_seg, freq_mhz, nominal_z, butt_x = baselines[a["id"]]
        dz = a["heightM"] - nominal_z
        # load_base() always recenters the boom to its geometric midpoint at
        # X=0 (the mast axis); for a butt-mounted instance, re-reference the
        # boom so its rear tip - near the reflector, see load_base()'s
        # comment - sits on the mast axis instead, before heading rotation
        # swings the offset boom around it. Height (Z) is unaffected either
        # way - every baseline's boom runs along X at constant Z.
        base_wires = shift_wires_x(wires, -butt_x) if a.get("mount") == "butt" else wires
        placed = gcp.transform_wires(base_wires, dz=dz, heading_deg=a["headingDeg"])
        participants.append({
            "id": a["id"], "model": a["model"], "polarization": a["polarization"],
            "headingDeg": a["headingDeg"], "wires": placed,
            "fedTag": offset + fed_tag, "fedSeg": fed_seg, "wireCount": len(placed),
        })
        offset += len(placed)

    joint_wires = list(itertools.chain(*[p["wires"] for p in participants]))

    deltas = {}
    for fed in participants:
        elev_deg, phi_deg, role_grid = solve_pattern(joint_wires, fed["fedTag"], fed["fedSeg"], ref_freq, fed["headingDeg"])
        iso = get_isolated(fed["model"], fed["polarization"])
        assert elev_deg == iso["elevDeg"] and phi_deg == iso["phiDeg"]
        nE, nP = len(elev_deg), len(phi_deg)
        delta = [[round(role_grid[ti][pi] - iso["gainDb"][ti][pi], 2) for pi in range(nP)] for ti in range(nE)]
        deltas[str(fed["id"])] = delta

    return {"ok": True, "freqMhz": ref_freq, "coarseElevDeg": elev_deg, "coarsePhiDeg": phi_deg, "deltas": deltas}


class Handler(http.server.SimpleHTTPRequestHandler):
    def do_GET(self):
        if self.path.rstrip("/") == "/compute-coupling":
            models = [[m, p] for (m, p) in BASELINE_FILES if os.path.isfile(BASELINE_FILES[(m, p)])]
            self._send_json(200, {"available": True, "models": models})
            return
        super().do_GET()

    def do_POST(self):
        if self.path.rstrip("/") != "/compute-coupling":
            self.send_error(404)
            return
        try:
            length = int(self.headers.get("Content-Length", 0))
            body = json.loads(self.rfile.read(length) or b"{}")
            antennas = body.get("antennas")
            if not isinstance(antennas, list):
                raise SolveError("bad_request", "Missing 'antennas' array")
            for a in antennas:
                for field in ("id", "model", "polarization", "heightM", "headingDeg"):
                    if field not in a:
                        raise SolveError("bad_request", "Antenna entry missing '{}'".format(field))
            result = compute_coupling(antennas)
            self._send_json(200, result)
        except SolveError as e:
            self._send_json(422, {"ok": False, "error": e.code, "message": e.message, "detail": e.detail})
        except (ValueError, json.JSONDecodeError) as e:
            self._send_json(400, {"ok": False, "error": "bad_request", "message": str(e)})
        except Exception as e:
            traceback.print_exc()
            self._send_json(500, {"ok": False, "error": "internal_error", "message": str(e)})

    def _send_json(self, status, payload):
        body = json.dumps(payload).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def log_message(self, fmt, *args):
        sys.stderr.write("%s - %s\n" % (self.address_string(), fmt % args))


def main():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 8765
    # Defaults to loopback-only for local dev; set COUPLING_SERVER_HOST=0.0.0.0
    # to accept connections from outside the host (e.g. behind a reverse proxy).
    host = os.environ.get("COUPLING_SERVER_HOST", "127.0.0.1")
    server = http.server.HTTPServer((host, port), Handler)
    print("Serving Detection Range Viewer + /compute-coupling on http://{}:{}/".format(host, port))
    print("Open http://{}:{}/index.html".format("127.0.0.1" if host == "0.0.0.0" else host, port))
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()

# Detection Range Viewer

A 3D viewer that turns NEC-2 antenna radiation patterns into estimated Motus
tag detection range, in kilometers, for several ground/terrain conditions.

## Usage

1. Run the parser against a folder of NEC-2 `.out` files (the `.nec` files
   don't need to be there, just their `.out` results):

   ```
   python parse_nec.py "C:\path\to\some\antenna\folder"
   ```

   This writes `patterns.json` next to `parse_nec.py` (pass a second argument
   to write it somewhere else).

2. Open `index.html` in a browser (or re-publish it as a Claude Artifact
   alongside the new `patterns.json`).

Re-run step 1 any time `.out` files in the source folder change, then reload
`index.html` — no code edits needed for a new antenna, as long as its `.out`
files follow the same conventions as PLC1669's (see below).

## What a `.out` file needs to be picked up

- A `RADIATION PATTERNS` table (NEC-2's standard `RP` output), swept over a
  full theta/phi grid (at least 5x5) — a single 1D E/H-plane cut is skipped.
- A `GN` ground-plane data card if it's a terrain case. A file with **no**
  `GN` card is instead treated as that antenna/polarization's **free-space
  baseline** (bucketed under the key `freeSpace`), which the viewer's live
  height/terrain model (see below) builds on.
- Optionally, a `GROUND TYPE: <label> (...)` comment card (`CM` line) in the
  source `.nec` file, to give the terrain a friendly name in the UI instead
  of falling back to its eps/sigma or filename.

The parser auto-detects which end of the theta sweep is straight up (zenith)
by finding the end where gain doesn't vary with phi — so it doesn't matter
whether a given NEC run swept theta as 0..90, -90..0, or something else. When
theta=0.0 is itself one of the swept samples, it's used directly as zenith
(NEC-2 always measures theta from +Z, and every model here is built with
Z=up), since the phi-invariance heuristic alone can be fooled by a real
on-axis pattern null — e.g. a vertical (Z-oriented) dipole element has a deep
null straight up along its own axis, where dB values swing sharply with phi
even though the underlying field is smooth.

## Generating `.out` files without 4nec2

`run_nec.py` runs a `.nec` card deck locally via [PyNEC](https://pypi.org/project/PyNEC/)
(`pip install PyNEC`) and writes a `.out` file in the same format `parse_nec.py`
expects, for the small set of cards this project's decks use (`CM`, `CE`,
`GW`, `GE`, `GN`, `EX`, `FR`, `RP`, `EN`):

```
python run_nec.py path\to\antenna.nec
```

Validated against this project's existing 4nec2-generated `.out` files (peak
gain matches within ~0.3 dB — normal numerical variation between NEC-2
implementations, not a modeling error). 4nec2 itself still works fine too, and
is the only option for cards `run_nec.py` doesn't handle.

PyNEC's `get_gain_tot()` returns its flat gain array ordered phi-major/
theta-minor (phi is the outer loop); `run_nec.py` reshapes it accordingly
(`reshape(n_phi, n_theta).T`). An earlier version reshaped it the other way
round, which doesn't error (same total element count) but silently produces
a geometrically scrambled pattern - own peak magnitude looks plausible, but
at the wrong theta/phi. If a future `.out` from this script has its peak in
an implausible direction, suspect this class of bug first.

For a new antenna model going forward, only **one** full-grid free-space run
per polarization is needed (via `run_nec.py` or 4nec2, `GE 0` and no `GN`
card) — height and terrain are handled live in the browser by an analytic
ground-reflection model built on that baseline, rather than needing a
separate NEC-2 run per height/terrain combination. Use a fine theta (elevation)
step on that RP card - 1 degree, matching the terrain runs - not the coarser
2.5 degree step this project used at first: the two-ray interference lobes
the live model produces get narrow at low grazing angles for tall antennas
or high frequencies, and a coarse baseline grid can miss the true peak by
more than a dB.

## Files

- `parse_nec.py` — the parser/generator. Self-contained, no dependencies
  beyond the Python standard library.
- `run_nec.py` — optional local NEC-2 runner (`.nec` -> `.out`) via PyNEC, for
  when 4nec2 isn't handy.
- `index.html` — the viewer. Reads `patterns.json` at load time.
- `patterns.json` — generated output, not hand-edited.

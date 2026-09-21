# Occult hybrid mode switches and barrier-augmented metastatic edge-rate identifiability

**Thesis #30.** Computational research, set out in Nile University B.Sc. chapter order for handoff.

**Depends on:** Thesis #20 (occult modes under partial observer) and Thesis #21 (metastasis graph barrier conductances).

**Author:** Kelechi Emeka Ogbonna  
**Email:** kelechiogbonna300@gmail.com  
**GitHub:** https://github.com/cloudynirvana  
**Date:** 21 September 2026

When occult hybrid mode switches and edge-wise desmoplastic conductances coexist on one anatomical graph toy, which metastatic edge rates stay practically identifiable from site-plus-barrier schedules that only see a mode-blind or mode-partial map?

On a three-node, two-edge toy they do not stay from site burdens alone. Structural rank is 4 of 6, and every shedding rate has an unbounded marginal variance. A duration-weighted stall on each edge, with no mode tag, restores local practical rank 6 of 6. Under that local rule every shedding rate stays: the widest marginal relative standard error is 0.1718, on the transport-limited active rate. The profile does not treat that rate as tightly as the local number does. On the mode-blind pool, and on a mode-partial map that resolves only the other edge, multipliers 0.78 through 1.40 of that rate remain inside the χ² cut. Recording that edge's stall in both modes shortens the inside set to 0.86 through 1.16. The shedding-limited active rate already has both neighbouring coarse-grid samples outside the cut on the blind pool (weighted residuals 32.16 and 20.92).

A continuous-flow alias fitted to the same blind record is practical rank 4 of 4 and lands on neither hybrid shedding rate (0.04544 against 0.0800 and 0.0280; 0.03699 against 0.0500 and 0.0250), at weighted cost 1826.

No number is taken from either parent deposit's results file. Thesis #20's pause contrast and Thesis #21's six-parameter ranks are not rows of this toy. This deposit does not claim a staging rule or a cure.

This is research only. It is not a medical device, not clinical decision support, not a dose, and not a cure. No document DOI is registered.

See [DISCLAIMER.md](DISCLAIMER.md). The manuscript is [THESIS.md](THESIS.md).

## Files

| Path | Role |
| --- | --- |
| `THESIS.md` | Manuscript (Chapters 1 to 5, Vancouver citations) |
| `THESIS.pdf` | PDF built from the Markdown |
| `build_pdf.py` | Regenerates `THESIS.pdf` |
| `CITATION.cff` | Citation metadata, no document DOI |
| `DISCLAIMER.md` | Research-only boundary |
| `sim/joint_toy.py` | Shared toy: hybrid switch, barrier schedules, Fisher ranks (seed 20260921) |
| `sim/results.json` | Numbers cited in Chapter Four |
| `sim/figures/` | Paths, spectra, marginal errors, profiles |

## Reproduce

```bash
python3 -m pip install -r sim/requirements.txt
python3 sim/joint_toy.py
python3 build_pdf.py
```

NumPy, SciPy and Matplotlib are required for the toy. The PDF step also needs the `markdown` and `weasyprint` packages. Regenerating the script rewrites `sim/results.json` and `sim/figures/`.

## Cite

Ogbonna KE. Occult hybrid mode switches and barrier-augmented metastatic edge-rate identifiability [Internet]. Thesis #30 computational research thesis. 21 September 2026 [cited YYYY Mon DD]. Available from: https://github.com/cloudynirvana/thesis-30-occult-modes-barrier-edge-rates

Machine-readable fields are in `CITATION.cff`. Add a document DOI there only after one exists.

Hub index, for cataloguing only: [research-theses-hub](https://github.com/cloudynirvana/research-theses-hub).

## Licence

Text and sketch code are MIT, with attribution. Computational research only.

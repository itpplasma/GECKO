# gecko.coils

Coil visualization and mesh generation.

```{eval-rst}
.. automodule:: gecko.coils
   :members:
   :undoc-members:
   :show-inheritance:
```

## STEP winding packs

Install the optional OpenCASCADE binding with `pip install -e '.[step]'`.
Import a STEP file and write MAKEGRID filaments with signed currents in amperes:

```sh
gecko coils-from-step winding_packs.step -o coils.cad \
    --currents 1200 -1200 --n-points 1024
```

One current broadcasts to all solids; otherwise supply one per solid in STEP
traversal order. Each solid becomes a separate group. `--reverse 0 2` reverses
the point ordering of solids 0 and 2. CAD has no physical current direction:
by default the largest component of each polygon's oriented area vector is
positive, and the loop starts at its lexicographically smallest section center.
Check this convention against the coil circuit before supplying signed currents.
No symmetry copies are created (`periods 1`, `mirror NIL`).

```python
from gecko.step import extract_step_centerlines, write_step_coils

lines_m = extract_step_centerlines("winding_packs.step", n_points=1024)
coils = write_step_coils("winding_packs.step", "coils.cad", currents=1200)
```

`--units mm` (default) and `--units m` select the CAD kernel's working unit.
Both respect the unit declared in the STEP header and return meters, following
[CadQuery's STEP unit conversion](https://cadquery.readthedocs.io/en/latest/importexport.html#importing-step).
This is not an override for incorrectly labeled STEP files. The [MAKEGRID output](https://princetonuniversity.github.io/STELLOPT/MAKEGRID.html)
contains `N` points plus an explicit repeated closing point with zero current and group
tag, followed by `end`. Zero-current coils remain in the Python result but are
omitted from the file; all-zero output is rejected. libneo treats zero-current
rows as separators.

The supported geometry is a closed chain of planar rectangular sections joined
by quadrilateral loft faces, such as ALPES winding packs. The importer separates
the two edge families of the boundary grid, recovers the shorter planar section
loops, takes their area centroids, and resamples the resulting closed polygon
uniformly in arc length. `n_points` counts distinct samples; the returned arrays
have shape `(n_points + 1, 3)` and repeat the first point. Increasing this count
reduces resampling error but does not reconstruct curvature absent from the CAD
section stations. Open shells, capped conductors, self-intersections and
unsupported topology raise `StepCoilError` with a coil index when available.
Smooth periodic sweeps (including B-spline faces) without section edges,
filleted sections, hollow packs and arbitrary solid medial axes are unsupported.
The regression suite explicitly checks that a valid closed rectangular B-spline
sweep raises `StepCoilError`. A section's area centroid defines the recovered
filament; CAD alone cannot identify a different intended filament for an
off-center conductor. Multi-filament cross-sections are not implemented.

The round-trip test uses the complete 160-point ALPES coil 0, its original
40 × 67 mm rectangular pack and an adapted loft/union construction (attribution
and licenses in [the fixture README](../../python/tests/data/alpes/README.md)).
It does not execute the upstream ALPES export script. Because the generated
sections are centered on the source vertices, section-center extraction should
recover those vertices within **1e-9 m**; this is measured before resampling.
Uniform arc-length resampling introduces new chords that can cut across the
original polygon's corners. The test bounds source-to-resampled distance,
sampled at 32 positions per source segment, by **0.2 mm**, and separately checks
that resampled-to-source distance stays below **1e-9 m**.

The test-only field oracle uses the exact finite straight-segment Biot–Savart
integral, with orientations matched. It reports extraction error (section-center
polygon versus source polygon) separately from resampling error (resampled
polygon versus section-center polygon). Relative vector errors must stay below
**1e-9** for extraction and **0.2%** for both resampling and the full round trip
at 27 points in `[-0.08, 0.08] m`, each more than 30 mm from the filament.
These fixture-specific checks do not validate a smooth physical coil field or
arbitrary CAD.

To save figure-ready comparisons and statistics, run:

```sh
GECKO_STEP_ARTIFACT_DIR=artifacts/step \
    pytest python/tests/test_step.py::test_alpes_roundtrip -v
```

This writes `centerline.csv` (resampled vertices paired with nearest source
polyline projections and forward distances, in meters), `field_errors.csv`
(positions, source/section-center/resampled fields in tesla, and separate
extraction/resampling/total relative errors), and `statistics.csv` (maximum
section-center extraction error, both directed resampling distances, separate
field errors and tolerances). The backward distance is sampled, not an exact
Hausdorff distance. STEP output remains in pytest's temporary directory.

```{eval-rst}
.. automodule:: gecko.step
   :members:
```

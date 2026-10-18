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
Smooth periodic sweeps without section edges, fillets, hollow packs and arbitrary
solid medial axes are currently unsupported. Multi-filament cross-sections are
not implemented.

The round-trip test uses the complete 160-point ALPES coil 0, its original
40 × 67 mm rectangular pack and loft/union construction (attribution and licenses
in `python/tests/data/README.md`). The acceptance thresholds are a **0.2 mm**
maximum bidirectional sampled distance between the original and extracted
polylines, and **0.2%** maximum relative Biot–Savart vector error at 27 points
in the cube `[-0.08, 0.08] m` (each more than 30 mm from the filament). Fields
use the exact finite straight-segment integral, with source and extracted
orientations matched. These tolerances apply to this fixture, not arbitrary CAD.

To save figure-ready comparisons and statistics, run:

```sh
GECKO_STEP_ARTIFACT_DIR=artifacts/step \
    pytest python/tests/test_step.py::test_alpes_roundtrip -v
```

This writes `centerline.csv` (extracted samples paired with nearest original
polyline projections, all in meters), `field_errors.csv` (positions, both fields
in tesla and relative errors), and `statistics.csv` (max deviation, max/RMS field
error and tolerances). STEP output remains in pytest's temporary directory.

```{eval-rst}
.. automodule:: gecko.step
   :members:
```

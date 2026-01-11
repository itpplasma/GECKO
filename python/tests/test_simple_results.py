from __future__ import annotations

import numpy as np

import pytest


def test_load_simple_wall_hits_units_cm_to_m(tmp_path) -> None:
    netCDF4 = pytest.importorskip("netCDF4")

    path = tmp_path / "results.nc"
    with netCDF4.Dataset(path, "w") as ds:
        ds.createDimension("particle", 4)
        ds.createDimension("xyz", 3)

        wall_hit = ds.createVariable("wall_hit", "i1", ("particle",))
        wall_hit[:] = np.array([0, 1, 0, 1], dtype=np.int8)

        wall_hit_cart = ds.createVariable("wall_hit_cart", "f8", ("xyz", "particle"))
        wall_hit_cart.units = "cm"
        wall_hit_cart[:, :] = np.array(
            [
                [0.0, 100.0, 0.0, 200.0],
                [0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0],
            ],
            dtype=float,
        )

        cos_inc = ds.createVariable("wall_hit_cos_incidence", "f8", ("particle",))
        cos_inc[:] = np.array([0.0, 0.1, 0.0, 0.2], dtype=float)

    from gecko.simple_results import load_simple_wall_hits

    hits = load_simple_wall_hits(path)
    assert hits.points_m.shape == (2, 3)
    assert np.allclose(hits.points_m[:, 0], np.array([1.0, 2.0]))
    assert hits.cos_incidence is not None
    assert np.allclose(hits.cos_incidence, np.array([0.1, 0.2]))

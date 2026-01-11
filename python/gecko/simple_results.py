from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class SimpleWallHits:
    points_m: np.ndarray  # (N,3)
    cos_incidence: np.ndarray | None  # (N,)
    angle_rad: np.ndarray | None  # (N,)


def load_simple_wall_hits(results_nc: Path | str) -> SimpleWallHits:
    """
    Load wall-hit points from SIMPLE `results.nc`.

    Conventions:
    - `wall_hit_cart` is stored in chartmap Cartesian units, written as a NetCDF
      `units` attribute (typically `cm` or `m`).
    - This loader returns points in **meters**.
    """
    from netCDF4 import Dataset

    results_nc = Path(results_nc)
    with Dataset(results_nc, "r") as ds:
        if "wall_hit" not in ds.variables or "wall_hit_cart" not in ds.variables:
            raise ValueError("results.nc missing wall_hit and/or wall_hit_cart")

        wall_hit = np.asarray(ds.variables["wall_hit"][:], dtype=np.int8)
        hit_cart = np.asarray(ds.variables["wall_hit_cart"][:], dtype=float)

        units = getattr(ds.variables["wall_hit_cart"], "units", "")
        if not units:
            units = getattr(ds, "cart_units", "")
        units = str(units).strip().lower()
        if units in ("", "cm"):
            cart_scale_to_m = 1.0e-2
        elif units == "m":
            cart_scale_to_m = 1.0
        else:
            raise ValueError(f"Unsupported wall_hit_cart units={units!r} (expected cm or m)")

        if hit_cart.ndim != 2:
            raise ValueError("wall_hit_cart must be 2D")
        # Accept both (3, N) and (N, 3) conventions.
        if hit_cart.shape[0] == 3:
            hit_cart = hit_cart
        elif hit_cart.shape[1] == 3:
            hit_cart = hit_cart.T
        else:
            raise ValueError("wall_hit_cart must have shape (3, N) or (N, 3)")

        mask = wall_hit != 0
        points_cart = hit_cart[:, mask].T
        points_m = cart_scale_to_m * points_cart

        cos_incidence = None
        angle_rad = None
        if "wall_hit_cos_incidence" in ds.variables:
            cos_incidence = np.asarray(ds.variables["wall_hit_cos_incidence"][:], dtype=float)[
                mask
            ]
        if "wall_hit_angle_rad" in ds.variables:
            angle_rad = np.asarray(ds.variables["wall_hit_angle_rad"][:], dtype=float)[mask]

    return SimpleWallHits(points_m=points_m, cos_incidence=cos_incidence, angle_rad=angle_rad)

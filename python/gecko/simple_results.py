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
    - `wall_hit_cart` is stored in **cm** in SIMPLE output.
    - This loader returns points in **meters**.
    """
    from netCDF4 import Dataset

    results_nc = Path(results_nc)
    with Dataset(results_nc, "r") as ds:
        if "wall_hit" not in ds.variables or "wall_hit_cart" not in ds.variables:
            raise ValueError("results.nc missing wall_hit and/or wall_hit_cart")

        wall_hit = np.asarray(ds.variables["wall_hit"][:], dtype=np.int8)
        hit_cart_cm = np.asarray(ds.variables["wall_hit_cart"][:], dtype=float)

        if hit_cart_cm.ndim != 2 or hit_cart_cm.shape[0] != 3:
            raise ValueError("wall_hit_cart must have shape (3, N)")

        mask = wall_hit != 0
        points_cm = hit_cart_cm[:, mask].T
        points_m = 0.01 * points_cm

        cos_incidence = None
        angle_rad = None
        if "wall_hit_cos_incidence" in ds.variables:
            cos_incidence = np.asarray(ds.variables["wall_hit_cos_incidence"][:], dtype=float)[
                mask
            ]
        if "wall_hit_angle_rad" in ds.variables:
            angle_rad = np.asarray(ds.variables["wall_hit_angle_rad"][:], dtype=float)[mask]

    return SimpleWallHits(points_m=points_m, cos_incidence=cos_incidence, angle_rad=angle_rad)


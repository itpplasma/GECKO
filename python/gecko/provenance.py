"""Version provenance for geometry exported from GECKO."""

from __future__ import annotations

import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import subprocess


def geometry_version() -> dict[str, str | None]:
    """Return the installed version and source revision, when available.

    Resolve Git relative to this package, never the caller's working directory.
    Wheels retain the revision captured by the package build hook.
    """
    root = Path(__file__).resolve().parents[2]
    try:
        from ._build_revision import GIT_DESCRIBE
    except ImportError:
        GIT_DESCRIBE = None
    describe = GIT_DESCRIBE
    if (root / "python" / "gecko" / "provenance.py").is_file() and (
        root / ".git"
    ).exists():
        try:
            describe = subprocess.check_output(
                ["git", "describe", "--tags", "--always", "--dirty"],
                cwd=root,
                stderr=subprocess.DEVNULL,
                text=True,
                timeout=5,
            ).strip()
        except (OSError, subprocess.SubprocessError):
            pass
    return {"gecko_version": version("gecko"), "git_describe": describe}


def write_stl_provenance(path: Path) -> None:
    """Write a version sidecar tied to the exported STL bytes."""
    provenance = geometry_version()
    provenance["stl_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    path.with_suffix(path.suffix + ".json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
    )

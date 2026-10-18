"""Retain the source revision in wheels built from a Git checkout."""

from pathlib import Path
import subprocess

from hatchling.builders.hooks.plugin.interface import BuildHookInterface


class CustomBuildHook(BuildHookInterface):
    def initialize(self, version, build_data):
        root = Path(self.root)
        revision = None
        if (root / ".git").exists():
            try:
                revision = subprocess.check_output(
                    ["git", "describe", "--tags", "--always", "--dirty"],
                    cwd=root,
                    text=True,
                    stderr=subprocess.DEVNULL,
                    timeout=5,
                ).strip()
            except (OSError, subprocess.SubprocessError):
                pass
        stamp = root / "build" / "gecko" / "_build_revision.py"
        stamp.parent.mkdir(parents=True, exist_ok=True)
        stamp.write_text(f"GIT_DESCRIBE = {revision!r}\n", encoding="utf-8")
        build_data["force_include"][str(stamp)] = "gecko/_build_revision.py"

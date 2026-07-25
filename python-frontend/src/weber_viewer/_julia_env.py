"""Bootstrap the scratch Julia project that :mod:`juliacall` runs against.

``juliacall`` resolves its Julia dependencies (PythonCall, OpenSSL_jll, ...) into
whatever project ``PYTHON_JULIAPKG_PROJECT`` points at. Pointing it at the
repository root would mutate the package's own ``Project.toml``, so we hand it a
dedicated scratch project that ``dev``s the local package instead.

This runs from Python rather than a wrapper script so ``uv run`` is the only
entry point users need.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

__all__ = ["ensure_julia_env", "julia_executable", "repo_root", "julia_env_dir"]

# src/weber_viewer/_julia_env.py -> src/weber_viewer -> src -> python-frontend
_FRONTEND_ROOT = Path(__file__).resolve().parents[2]


def repo_root() -> Path:
    """Absolute path to the WeberElectrodynamics.jl repository root."""
    return _FRONTEND_ROOT.parent


def julia_env_dir() -> Path:
    """Scratch Julia project directory (gitignored, regenerated on demand)."""
    override = os.environ.get("WEBER_VIEWER_JULIA_ENV")
    if override:
        return Path(override).expanduser().resolve()
    return _FRONTEND_ROOT / ".julia-env"


def julia_executable() -> str:
    """Locate a real ``julia`` binary, resolving the juliaup launcher shim.

    ``juliacall`` invokes the executable directly and the juliaup launcher does
    not always cooperate, so resolve the shim to the concrete install backing
    juliaup's *default* channel. Asking Julia for its own ``Sys.BINDIR`` beats
    globbing the depot, which would happily pick a newer prerelease over the
    release channel the user actually selected.
    """
    explicit = os.environ.get("PYTHON_JULIAPKG_EXE")
    if explicit:
        return explicit

    julia = shutil.which("julia")
    if julia is None:
        raise RuntimeError(
            "'julia' not found on PATH — install juliaup or add julia to PATH"
        )

    if not os.path.realpath(julia).endswith("julialauncher"):
        return julia

    try:
        result = subprocess.run(
            [julia, "--startup-file=no", "-e", "print(Sys.BINDIR)"],
            capture_output=True,
            text=True,
            check=True,
            timeout=120,
        )
    except (subprocess.SubprocessError, OSError):
        return julia

    resolved = Path(result.stdout.strip()) / "julia"
    if resolved.is_file() and os.access(resolved, os.X_OK):
        return str(resolved)
    return julia


def ensure_julia_env() -> Path:
    """Create the scratch Julia project if absent and export juliacall's env vars.

    Must be called before ``juliacall`` is imported. Returns the project path.
    """
    env_dir = julia_env_dir()
    julia_bin = julia_executable()

    if not (env_dir / "Project.toml").exists():
        env_dir.mkdir(parents=True, exist_ok=True)
        print(f"[weber-viewer] creating Julia project at {env_dir}")
        subprocess.run(
            [
                julia_bin,
                f"--project={env_dir}",
                "-e",
                f'using Pkg; Pkg.develop(path=raw"{repo_root()}")',
            ],
            check=True,
        )

    os.environ.setdefault("PYTHON_JULIAPKG_PROJECT", str(env_dir))
    os.environ.setdefault("PYTHON_JULIAPKG_EXE", julia_bin)
    return env_dir

#!/usr/bin/env bash
# Run a weber-viewer example.
#
# Usage:
#   ./python-frontend/run-example.sh                         # two_body_streaming.py (default)
#   ./python-frontend/run-example.sh path/to/other.py        # any other script
#   ./python-frontend/run-example.sh examples/other.py
#
# First run bootstraps two isolated environments:
#   .venv-viewer/            — Python venv with weber-viewer installed (editable)
#   .venv-viewer/julia-env/  — dedicated Julia project that `dev`s the local
#                              WeberElectrodynamics package, so juliacall can
#                              freely add PythonCall/OpenSSL_jll here without
#                              mutating the package's own Project.toml.

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." &>/dev/null && pwd)"
VENV="${REPO_ROOT}/.venv-viewer"
JULIA_ENV="${VENV}/julia-env"
DEFAULT_EXAMPLE="${SCRIPT_DIR}/examples/two_body_streaming.py"

# ---- Locate julia -----------------------------------------------------------
if ! JULIA_BIN="$(command -v julia 2>/dev/null)"; then
    echo "error: 'julia' not found on PATH — install juliaup or add julia to PATH" >&2
    exit 1
fi

# ---- Bootstrap Python venv + install once ----------------------------------
if [[ ! -x "${VENV}/bin/python" ]]; then
    echo "[run-example] creating Python venv at ${VENV}"
    python3 -m venv "${VENV}"
    "${VENV}/bin/pip" install --quiet --upgrade pip
fi

if ! "${VENV}/bin/python" -c "import weber_viewer" 2>/dev/null; then
    echo "[run-example] installing weber-viewer (editable)"
    "${VENV}/bin/pip" install --quiet -e "${SCRIPT_DIR}"
fi

# ---- Bootstrap dedicated Julia project -------------------------------------
# Keeps juliapkg's PythonCall/OpenSSL_jll additions out of the package's own
# Project.toml by giving juliacall a separate scratch project that dev's the
# local package.
if [[ ! -f "${JULIA_ENV}/Project.toml" ]]; then
    echo "[run-example] creating Julia project at ${JULIA_ENV}"
    mkdir -p "${JULIA_ENV}"
    "${JULIA_BIN}" --project="${JULIA_ENV}" -e "
        using Pkg
        Pkg.develop(path=raw\"${REPO_ROOT}\")
    "
fi

# ---- Point juliacall at the dedicated project ------------------------------
export PYTHON_JULIAPKG_PROJECT="${JULIA_ENV}"
export PYTHON_JULIAPKG_EXE="${JULIA_BIN}"

# ---- Resolve script to run --------------------------------------------------
SCRIPT="${1:-${DEFAULT_EXAMPLE}}"
shift || true

if [[ ! -f "${SCRIPT}" && "${SCRIPT}" != -* ]]; then
    for CANDIDATE in \
        "${SCRIPT_DIR}/${SCRIPT}" \
        "${SCRIPT_DIR}/examples/${SCRIPT}" \
        "${SCRIPT_DIR}/examples/$(basename "${SCRIPT}")"; do
        if [[ -f "${CANDIDATE}" ]]; then
            SCRIPT="${CANDIDATE}"
            break
        fi
    done
fi

echo "[run-example] running ${SCRIPT}"
exec "${VENV}/bin/python" "${SCRIPT}" "$@"

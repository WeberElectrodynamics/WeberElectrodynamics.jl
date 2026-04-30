#!/usr/bin/env bash
# Run a weber-viewer example.
#
# Usage:
#   ./python-frontend/run-example.sh                         # two-body (default)
#   ./python-frontend/run-example.sh --list                  # list bundled examples
#   ./python-frontend/run-example.sh smoke --steps 50        # headless bridge probe
#   ./python-frontend/run-example.sh three-body              # 2D three-body showcase
#   ./python-frontend/run-example.sh path/to/other.py        # any other script
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
THREE_BODY_EXAMPLE="${SCRIPT_DIR}/examples/three_body_polygon.py"
SMOKE_EXAMPLE="${SCRIPT_DIR}/examples/smoke_probe.py"

usage() {
    cat <<EOF
Usage:
  ./python-frontend/run-example.sh [example] [args...]

Examples:
  two-body    3D two-body streaming viewer (default)
  three-body  2D three-body polygon viewer
  smoke       Headless Julia bridge probe

You may also pass a Python script path.
EOF
}

list_examples() {
    cat <<EOF
two-body    ${DEFAULT_EXAMPLE}
three-body  ${THREE_BODY_EXAMPLE}
smoke       ${SMOKE_EXAMPLE}
EOF
}

case "${1:-}" in
    --help|-h)
        usage
        exit 0
        ;;
    --list)
        list_examples
        exit 0
        ;;
esac

# ---- Resolve script to run --------------------------------------------------
SCRIPT_ARG="${1:-two-body}"
if [[ $# -gt 0 ]]; then
    shift
fi

case "${SCRIPT_ARG}" in
    two-body|two_body|two_body_streaming|two_body_streaming.py)
        SCRIPT="${DEFAULT_EXAMPLE}"
        ;;
    three-body|three_body|three_body_polygon|three_body_polygon.py)
        SCRIPT="${THREE_BODY_EXAMPLE}"
        ;;
    smoke|smoke_probe|smoke_probe.py)
        SCRIPT="${SMOKE_EXAMPLE}"
        ;;
    *)
        SCRIPT="${SCRIPT_ARG}"
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
        if [[ ! -f "${SCRIPT}" || "${SCRIPT}" == -* ]]; then
            echo "error: unknown example or script: ${SCRIPT_ARG}" >&2
            echo >&2
            usage >&2
            echo >&2
            echo "Bundled examples:" >&2
            list_examples >&2
            exit 2
        fi
        ;;
esac

# ---- Locate julia -----------------------------------------------------------
if ! JULIA_BIN="$(command -v julia 2>/dev/null)"; then
    echo "error: 'julia' not found on PATH — install juliaup or add julia to PATH" >&2
    exit 1
fi
JULIA_TARGET="$(readlink "${JULIA_BIN}" 2>/dev/null || true)"
if [[ "${JULIA_TARGET}" == *julialauncher ]]; then
    JULIAUP_ROOT="${JULIAUP_DEPOT_PATH:-${HOME}/.julia}/juliaup"
    for CANDIDATE in \
        "${JULIAUP_ROOT}"/julia-[0-9]*.[0-9]*.[0-9]+*/bin/julia \
        "${JULIAUP_ROOT}"/julia-[0-9]*.[0-9]*.[0-9]+*/Julia-*.app/Contents/Resources/julia/bin/julia \
        "${JULIAUP_ROOT}"/julia-*/bin/julia \
        "${JULIAUP_ROOT}"/julia-*/Julia-*.app/Contents/Resources/julia/bin/julia; do
        if [[ -x "${CANDIDATE}" ]]; then
            JULIA_BIN="${CANDIDATE}"
            break
        fi
    done
fi

# ---- Bootstrap Python venv + install once ----------------------------------
if [[ ! -x "${VENV}/bin/python" ]]; then
    echo "[run-example] creating Python venv at ${VENV}"
    python3 -m venv "${VENV}"
    "${VENV}/bin/pip" install --quiet --upgrade pip
fi

if ! "${VENV}/bin/python" -c "import weber_viewer" 2>/dev/null || \
   [[ ! -x "${VENV}/bin/weber-viewer-smoke" ]]; then
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

echo "[run-example] running ${SCRIPT}"
"${VENV}/bin/python" "${SCRIPT}" "$@"

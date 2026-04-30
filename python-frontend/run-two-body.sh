#!/usr/bin/env bash
# Run the 3D two-body streaming viewer demo.

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
"${SCRIPT_DIR}/run-example.sh" two-body "$@"

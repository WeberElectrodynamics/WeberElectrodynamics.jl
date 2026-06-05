#!/usr/bin/env bash
# Run the 2D three-body polygon viewer demo.

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
"${SCRIPT_DIR}/run-example.sh" three-body "$@"

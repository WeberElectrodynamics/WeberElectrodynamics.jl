#!/usr/bin/env bash
# Run the headless Julia bridge smoke probe.

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
"${SCRIPT_DIR}/run-example.sh" smoke "$@"

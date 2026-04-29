#!/usr/bin/env bash
# Usage: ./release.sh [patch|minor|major] [--execute]
#
# Default: DRY RUN — shows what would happen, makes no changes.
# Pass --execute to actually bump, commit, push and trigger JuliaRegistrator.
#
# Before running, populate RELEASENOTES.md with the notes for this release.
# On --execute the notes are moved into CHANGELOG.md and RELEASENOTES.md is reset.
set -euo pipefail

# ── Parse arguments ────────────────────────────────────────────────────────────
BUMP=""
EXECUTE=false

for arg in "$@"; do
  case $arg in
    patch|minor|major) BUMP="$arg" ;;
    --execute) EXECUTE=true ;;
    *) echo "Usage: $0 [patch|minor|major] [--execute]"; exit 1 ;;
  esac
done

BUMP="${BUMP:-patch}"
REPO="WeberElectrodynamics/WeberElectrodynamics.jl"

# ── Read current version ───────────────────────────────────────────────────────
CURRENT=$(grep '^version' Project.toml | sed 's/version = "\(.*\)"/\1/')
IFS='.' read -r MAJOR MINOR PATCH <<< "$CURRENT"

case $BUMP in
  patch) PATCH=$((PATCH + 1)) ;;
  minor) MINOR=$((MINOR + 1)); PATCH=0 ;;
  major) MAJOR=$((MAJOR + 1)); MINOR=0; PATCH=0 ;;
esac

NEW="$MAJOR.$MINOR.$PATCH"
DATE=$(date +%Y-%m-%d)

# ── Validate RELEASENOTES.md ───────────────────────────────────────────────────
if [ ! -f RELEASENOTES.md ]; then
  echo "ERROR: RELEASENOTES.md not found."
  echo "Create it and add release notes before running this script."
  exit 1
fi

# Strip blank lines and HTML comments, including multi-line comment blocks;
# check something remains.
NOTES=$(python3 - <<'PYEOF'
import re
from pathlib import Path

content = Path("RELEASENOTES.md").read_text()
content = re.sub(r"<!--.*?-->", "", content, flags=re.S)
lines = [line.rstrip() for line in content.splitlines()]
while lines and not lines[0].strip():
    lines.pop(0)
while lines and not lines[-1].strip():
    lines.pop()
print("\n".join(lines))
PYEOF
)

if [ -z "$NOTES" ]; then
  echo "ERROR: RELEASENOTES.md is empty."
  echo "Add release notes to RELEASENOTES.md before releasing."
  exit 1
fi

NOTE_LINES=$(printf '%s' "$NOTES" | wc -l | tr -d ' ')

# ── Dry-run output ─────────────────────────────────────────────────────────────
if [ "$EXECUTE" = false ]; then
  echo ""
  echo "DRY RUN — pass --execute to apply"
  echo ""
  echo "  Current version : $CURRENT"
  echo "  New version     : $NEW"
  echo "  Release notes   : ($NOTE_LINES lines from RELEASENOTES.md)"
  echo ""
  echo "  Would modify    : Project.toml, CHANGELOG.md, RELEASENOTES.md"
  echo "  Would commit    : \"release: v$NEW\""
  echo "  Would push      : origin main"
  echo "  Would comment   : @JuliaRegistrator register on HEAD"
  echo ""
  echo "Release notes preview:"
  echo "────────────────────────────────────────"
  printf '%s\n' "$NOTES"
  echo "────────────────────────────────────────"
  echo ""
  exit 0
fi

# ── Execute ────────────────────────────────────────────────────────────────────
echo "Releasing $CURRENT → $NEW"

ACTIVE_GH_USER=$(gh auth status --hostname github.com --json hosts --jq '.hosts["github.com"][] | select(.active == true and .state == "success") | .login' 2>/dev/null || true)
if [ "$ACTIVE_GH_USER" != "WeberElectrodynamics" ]; then
  echo "ERROR: active gh account must be WeberElectrodynamics for JuliaRegistrator."
  echo "Current active account: ${ACTIVE_GH_USER:-none}"
  echo "Run: gh auth switch --user WeberElectrodynamics"
  exit 1
fi

# Full notes content (preserving blank lines between sections)
FULL_NOTES="$NOTES"

# ── Bump Project.toml ─────────────────────────────────────────────────────────
python3 - <<PYEOF
from pathlib import Path

path = Path("Project.toml")
content = path.read_text()
old = 'version = "$CURRENT"'
new = 'version = "$NEW"'
if old not in content:
    raise SystemExit("ERROR: Project.toml version bump failed; did not find " + old)
path.write_text(content.replace(old, new, 1))
PYEOF

# ── Update CHANGELOG.md ───────────────────────────────────────────────────────
# Inserts "## [X.Y.Z] - DATE\n\n<notes>" after the [Unreleased] heading
export FULL_NOTES NEW DATE
python3 - <<'PYEOF'
import os
from pathlib import Path

content = Path("CHANGELOG.md").read_text()
new_entry = f"## [Unreleased]\n\n## [{os.environ['NEW']}] - {os.environ['DATE']}\n\n{os.environ['FULL_NOTES']}"
new = content.replace("## [Unreleased]", new_entry, 1)
Path("CHANGELOG.md").write_text(new)
PYEOF

if ! grep -q "## \[$NEW\]" CHANGELOG.md; then
  echo "ERROR: CHANGELOG.md update failed"
  exit 1
fi
echo "Updated CHANGELOG.md"

# ── Reset RELEASENOTES.md ─────────────────────────────────────────────────────
cat > RELEASENOTES.md << 'TEMPLATE'
<!-- Add release notes here before running ./release.sh -->
<!-- This file is cleared automatically on each release execute -->
TEMPLATE
echo "Reset RELEASENOTES.md"

# ── Commit & push ─────────────────────────────────────────────────────────────
git add Project.toml CHANGELOG.md RELEASENOTES.md
git commit -m "release: v$NEW"
git push origin main

# ── Trigger JuliaRegistrator ──────────────────────────────────────────────────
SHA=$(git rev-parse HEAD)
BODY=$(printf '@JuliaRegistrator register\n\nRelease notes:\n\nSee CHANGELOG.md for full details.\n\n%s' "$FULL_NOTES")
gh api "/repos/$REPO/commits/$SHA/comments" -f body="$BODY"

echo ""
echo "✓ v$NEW pushed and JuliaRegistrator triggered"
echo "  Registry PR will merge in ~15 min"
echo "  TagBot will tag automatically after that"
echo "  Docs /stable will deploy on the new tag"
echo "  Papers.yml will compile PDFs once the GitHub Release is published"

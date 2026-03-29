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

# Strip blank lines and HTML comment lines; check something remains
NOTES=$(grep -v '^\s*$' RELEASENOTES.md | grep -v '^\s*<!--' || true)

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

# Full notes content (preserving blank lines between sections)
FULL_NOTES=$(python3 -c "
import re
content = open('RELEASENOTES.md').read()
# Strip leading/trailing blank lines and comment lines
lines = content.splitlines()
filtered = [l for l in lines if not l.strip().startswith('<!--')]
# Strip leading/trailing blank lines
while filtered and not filtered[0].strip():
    filtered.pop(0)
while filtered and not filtered[-1].strip():
    filtered.pop()
print('\n'.join(filtered))
")

# ── Bump Project.toml ─────────────────────────────────────────────────────────
sed -i '' "s/^version = \"$CURRENT\"/version = \"$NEW\"/" Project.toml

# ── Update CHANGELOG.md ───────────────────────────────────────────────────────
# Inserts "## [X.Y.Z] - DATE\n\n<notes>" after the [Unreleased] heading
python3 - <<PYEOF
import sys
notes = open('RELEASENOTES.md').read()
lines = notes.splitlines()
filtered = [l for l in lines if not l.strip().startswith('<!--')]
while filtered and not filtered[0].strip():
    filtered.pop(0)
while filtered and not filtered[-1].strip():
    filtered.pop()
clean_notes = '\n'.join(filtered)

content = open('CHANGELOG.md').read()
new_entry = f"## [Unreleased]\n\n## [$NEW] - $DATE\n\n{clean_notes}"
new = content.replace("## [Unreleased]", new_entry, 1)
open('CHANGELOG.md', 'w').write(new)
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

#!/usr/bin/env bash
# Usage: ./release.sh [patch|minor|major]
# - Bumps version in Project.toml
# - Moves [Unreleased] changelog entries to the new version
# - Commits, pushes, and triggers JuliaRegistrator automatically
set -euo pipefail

BUMP=${1:-patch}
REPO="WeberElectrodynamics/WeberElectrodynamics.jl"

# ── Read current version ───────────────────────────────────────────────────────
CURRENT=$(grep '^version' Project.toml | sed 's/version = "\(.*\)"/\1/')
IFS='.' read -r MAJOR MINOR PATCH <<< "$CURRENT"

case $BUMP in
  patch) PATCH=$((PATCH + 1)) ;;
  minor) MINOR=$((MINOR + 1)); PATCH=0 ;;
  major) MAJOR=$((MAJOR + 1)); MINOR=0; PATCH=0 ;;
  *) echo "Usage: $0 [patch|minor|major]"; exit 1 ;;
esac

NEW="$MAJOR.$MINOR.$PATCH"
DATE=$(date +%Y-%m-%d)

echo "Releasing $CURRENT → $NEW"

# ── Bump Project.toml ─────────────────────────────────────────────────────────
sed -i '' "s/^version = \"$CURRENT\"/version = \"$NEW\"/" Project.toml

# ── Update CHANGELOG.md ───────────────────────────────────────────────────────
# Inserts "## [X.Y.Z] - DATE" after the [Unreleased] line
python3 - <<EOF
content = open("CHANGELOG.md").read()
new = content.replace(
    "## [Unreleased]",
    "## [Unreleased]\n\n## [$NEW] - $DATE",
    1
)
open("CHANGELOG.md", "w").write(new)
EOF

if ! grep -q "## \[$NEW\]" CHANGELOG.md; then
  echo "ERROR: CHANGELOG.md update failed — [Unreleased] section not found"
  exit 1
fi
echo "Updated Project.toml and CHANGELOG.md"
grep -A1 "## \[$NEW\]" CHANGELOG.md

# ── Extract release notes (guard against empty section) ───────────────────────
NOTES=$(python3 -c "
import re
content = open('CHANGELOG.md').read()
m = re.search(r'## \[$NEW\][^\n]*\n(.*?)(?=\n## \[|\Z)', content, re.DOTALL)
print(m.group(1).strip() if m else '')
")

if [ -z "$NOTES" ]; then
  echo "ERROR: No release notes found under '## [$NEW]' in CHANGELOG.md"
  echo "Add content under the new version heading before releasing."
  exit 1
fi
echo "Release notes extracted ($(printf '%s' "$NOTES" | wc -l | tr -d ' ') lines)"

# ── Commit & push ─────────────────────────────────────────────────────────────
git add Project.toml CHANGELOG.md
git commit -m "release: v$NEW"
git push origin main

# ── Trigger JuliaRegistrator ──────────────────────────────────────────────────
SHA=$(git rev-parse HEAD)
BODY=$(printf '@JuliaRegistrator register\n\nRelease notes:\n\nSee CHANGELOG.md for full details.\n\n%s' "$NOTES")
gh api "/repos/$REPO/commits/$SHA/comments" -f body="$BODY"

echo ""
echo "✓ v$NEW pushed and JuliaRegistrator triggered"
echo "  Registry PR will merge in ~15 min"
echo "  TagBot will tag automatically after that"
echo "  Docs /stable will deploy on the new tag"

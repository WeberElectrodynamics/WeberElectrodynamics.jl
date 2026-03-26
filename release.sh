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

echo "Updated Project.toml and CHANGELOG.md"
grep -A1 "## \[$NEW\]" CHANGELOG.md

# ── Commit & push ─────────────────────────────────────────────────────────────
git add Project.toml CHANGELOG.md
git commit -m "release: v$NEW"
git push origin main

# ── Trigger JuliaRegistrator ──────────────────────────────────────────────────
SHA=$(git rev-parse HEAD)
gh api "/repos/$REPO/commits/$SHA/comments" -f body='@JuliaRegistrator register'

echo ""
echo "✓ v$NEW pushed and JuliaRegistrator triggered"
echo "  Registry PR will merge in ~15 min"
echo "  TagBot will tag automatically after that"
echo "  Docs /stable will deploy on the new tag"

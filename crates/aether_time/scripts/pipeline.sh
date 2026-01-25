#!/usr/bin/env bash
set -euo pipefail

# Simple local pipeline for Rust workspaces/crates.
# - Runs tests for all packages
# - Runs clippy (non-destructive; no --fix)
# - Builds docs for the workspace
#
# Usage:
#   bash scripts/pipeline.sh [--all-features] [--release] [--no-doc]
#
# Notes:
# - Defaults to debug builds without all features to be faster.
# - Set --all-features to expand coverage when desired.

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
root="$(cd "${here}/.." && pwd)"
cd "${root}"

ALL_FEATURES=false
RELEASE=false
BUILD_DOCS=true

while [[ $# -gt 0 ]]; do
  case "$1" in
    --all-features) ALL_FEATURES=true; shift ;;
    --release) RELEASE=true; shift ;;
    --no-doc) BUILD_DOCS=false; shift ;;
    *) echo "Unknown option: $1" >&2; exit 2 ;;
  esac
done

features=( )
${ALL_FEATURES} && features+=("--all-features")

profile=( )
${RELEASE} && profile+=("--release")

echo "==> cargo test ${features[*]} ${profile[*]} --workspace --all-targets"
cargo test --workspace --all-targets ${features[@]:-} ${profile[@]:-}

echo "==> cargo clippy ${features[*]} ${profile[*]} --workspace --all-targets -- -D warnings"
# Non-destructive: do NOT pass --fix
if ! cargo clippy --workspace --all-targets ${features[@]:-} ${profile[@]:-} -- -D warnings; then
  echo "clippy reported issues. Fix them and re-run." >&2
  exit 1
fi

if ${BUILD_DOCS}; then
  echo "==> cargo doc ${features[*]} --workspace --no-deps"
  cargo doc --workspace --no-deps ${features[@]:-}
fi

echo "Done."


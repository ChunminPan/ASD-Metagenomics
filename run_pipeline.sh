#!/usr/bin/env bash
# New publication entry point; see docs/PROVENANCE.md for historical sources.
set -Eeuo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
exec python3 "$repo_dir/scripts/pipeline.py" "$@"

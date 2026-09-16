#!/bin/sh
set -eu

ROOT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
MODEL_DIR="$ROOT_DIR/applications/regularisedHerschelBulkley"

if ! command -v wmake >/dev/null 2>&1; then
    echo "wmake not found. Source OpenFOAM Foundation v14 before building the custom viscosity model." >&2
    exit 1
fi

cd "$MODEL_DIR"
wmake libso

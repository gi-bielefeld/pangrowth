#!/usr/bin/env bash
set -euo pipefail

GGCAT_TAG="${GGCAT_TAG:-v2.0.0}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${GGCAT_BUILD_DIR:-/tmp/ggcat-${GGCAT_TAG}}"
OUT_DIR="$ROOT/third_party/ggcat-cpp-api/lib"

if [[ ! -d "$BUILD_DIR/.git" ]]; then
  git clone --depth 1 --branch "$GGCAT_TAG" https://github.com/algbio/ggcat.git "$BUILD_DIR"
fi

cargo build --release --manifest-path "$BUILD_DIR/Cargo.toml" --package ggcat-cpp-bindings
mkdir -p "$OUT_DIR"
cp "$BUILD_DIR/target/release/libggcat_cpp_bindings.a" "$OUT_DIR/"

interop="$(find "$BUILD_DIR/target/release/build" -name libggcat_cxx_interop.a -print -quit)"
if [[ -z "$interop" ]]; then
  echo "Could not find libggcat_cxx_interop.a under $BUILD_DIR/target/release/build" >&2
  exit 1
fi
cp "$interop" "$OUT_DIR/"

echo "Wrote ggcat static libraries to $OUT_DIR"

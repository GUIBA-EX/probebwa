#!/usr/bin/env bash
set -euo pipefail

# Upstream's .cargo/config.toml defaults to `-C target-cpu=native`, which is
# the right choice when compiling on the machine that will run the binary
# (this project's normal use case) but wrong here: bioconda builds once on
# its own build machine and ships the resulting binary to users on
# arbitrary, possibly older or different, CPUs. A target-cpu=native binary
# built here would crash with an illegal-instruction fault on any machine
# missing a feature the build machine happens to have. RUSTFLAGS overrides
# the config value entirely (it doesn't merge with it), so this pins the
# build to a portable baseline instead: x86-64-v2 (SSE4.2/POPCNT, universal
# on any x86_64 CPU from the last ~15 years) on x86_64, and the plain
# architecture default (no target-cpu override at all -- there's no
# equally-universal named "aarch64-v2" baseline) elsewhere.
case "$(uname -m)" in
    x86_64) export RUSTFLAGS="-C target-cpu=x86-64-v2" ;;
    *) export RUSTFLAGS="" ;;
esac

cargo install --locked --root "${PREFIX}" --path .

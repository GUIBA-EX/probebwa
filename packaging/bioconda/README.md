# Bioconda recipe (draft)

`meta.yaml`/`build.sh` here are a ready-to-submit Bioconda recipe for the
`v0.1.0` tag, verified locally by extracting the tag's GitHub-generated
source tarball and running `build.sh` against it end-to-end (`cargo install
--locked` succeeds, the installed binary runs `--help` correctly).

`meta.yaml`'s `source.sha256` is already the real sha256 of
`https://github.com/GUIBA-EX/probebwa/archive/refs/tags/v0.1.0.tar.gz` --
recompute it if you retarget this recipe at a different tag.

## To actually submit this to Bioconda

This part is a public PR to someone else's repository and isn't done here:

1. Fork `bioconda/bioconda-recipes` on GitHub.
2. Copy `meta.yaml` and `build.sh` from this directory into
   `recipes/probebwa/` in your fork (paths matter: Bioconda's tooling looks
   for `recipes/<name>/meta.yaml`).
3. Open a PR against `bioconda/bioconda-recipes`. Their CI lints and
   test-builds the recipe automatically; a Bioconda maintainer reviews and
   merges it, after which the package is built and uploaded to the
   `bioconda` conda channel automatically.

## Why `build.sh` overrides `RUSTFLAGS`

The main repo's `.cargo/config.toml` defaults to `-C target-cpu=native`,
which is correct for a normal from-source build (compile on the machine
that will run it) but wrong for a distributed binary: Bioconda builds once
on its own build machine and ships that exact binary to users on arbitrary
CPUs, so a `target-cpu=native` binary risks an illegal-instruction crash on
any machine missing a feature the build machine happened to have.
`build.sh` overrides it to a portable baseline instead (`x86-64-v2` on
x86_64; no override at all on other architectures, since there's no
equally-universal named baseline there).

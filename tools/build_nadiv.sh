#!/usr/bin/env bash
# Install nadiv (GPL >= 2; called as a library, never vendored) at a pinned
# CRAN version into .rlib/. Skips the install when .rlib/ already holds it.
set -euo pipefail
VERSION=2.18.0
SHA256=e6748d3527b8682c68d74e2657a696361aae12a1e6b681ed4c87b05891e31059
LIB="$PIXI_PROJECT_ROOT/.rlib"
if [[ "$(Rscript -e "cat(tryCatch(as.character(packageVersion('nadiv', lib.loc = '$LIB')), error = function(e) ''))")" == "$VERSION" ]]; then
  exit 0
fi
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
tarball="$tmp/nadiv_$VERSION.tar.gz"
# CRAN moves superseded versions to Archive/.
curl -sfL -o "$tarball" "https://cran.r-project.org/src/contrib/nadiv_$VERSION.tar.gz" ||
  curl -sfL -o "$tarball" "https://cran.r-project.org/src/contrib/Archive/nadiv/nadiv_$VERSION.tar.gz"
echo "$SHA256  $tarball" | sha256sum --check --quiet || {
  echo "nadiv $VERSION tarball does not match the pinned sha256; refusing to install" >&2
  exit 1
}
rm -rf "$LIB"
mkdir -p "$LIB"
R CMD INSTALL --no-docs --no-html --library="$LIB" "$tarball" >"$tmp/install.log" 2>&1 || {
  cat "$tmp/install.log" >&2
  exit 1
}

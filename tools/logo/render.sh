#!/bin/sh
# Rebuild man/figures/logo.svg (text outlined, so it needs no fonts) and logo.png.
# Needs python3, Inkscape >= 1.0, and the Comfortaa font (Bold).
set -e
here="$(cd "$(dirname "$0")" && pwd)"
figs="$here/../../man/figures"
tmp="$(mktemp --suffix=.svg)"
trap 'rm -f "$tmp"' EXIT
python3 "$here/make_logo.py" "$tmp"
inkscape "$tmp" --export-text-to-path --export-plain-svg --export-filename="$figs/logo.svg" 2>/dev/null
# 240 px wide (about 278 tall): the usethis::use_logo() size, sharp at the README's height="139".
inkscape "$figs/logo.svg" --export-type=png --export-width=240 --export-filename="$figs/logo.png" 2>/dev/null

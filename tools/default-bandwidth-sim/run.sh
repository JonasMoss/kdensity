#!/usr/bin/env bash
# Run every unfinished cell of sim.R with WORKERS parallel single-threaded
# workers under nice. Finished cells (results/<cell>.csv) are skipped, so the
# script can simply be restarted after an interruption.
#
#   KDENSITY_LIB=~/kdensity-sim/lib ./run.sh 16
set -euo pipefail
workers="${1:-4}"
here="$(cd "$(dirname "$0")" && pwd)"
outdir="${OUTDIR:-$here/results}"
mkdir -p "$outdir" "$here/logs"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
# Largest cells first so the tail of the run is short.
Rscript "$here/sim.R" list |
  awk -F'__' '{ n = substr($4, 2); print n "\t" $0 }' |
  sort -k1,1nr | cut -f2 |
  while read -r cell; do
    [ -f "$outdir/$cell.csv" ] || echo "$cell"
  done |
  xargs -P "$workers" -I{} sh -c \
    'nice -n 10 Rscript "$0/sim.R" run "$1" "$2" > "$0/logs/$1.log" 2>&1 || echo "FAILED $1"' \
    "$here" {} "$outdir"
echo "finished: $(ls "$outdir"/*.csv 2>/dev/null | wc -l) of $(Rscript "$here/sim.R" list | wc -l) cells"

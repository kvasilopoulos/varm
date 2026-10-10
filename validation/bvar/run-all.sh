#!/bin/sh
# Re-run all BVAR_ cross-validation scripts. Usage: VARM_SRC=<scratch copy of the repo> sh validation/bvar/run-all.sh
# (VARM_SRC avoids compiling inside the shared working tree; default is the repo itself.)
cd "$(dirname "$0")" || exit 1
for f in 01_minnesota 02_dynamics 03_proxy 04_signs 05_zerosign 06_hetero 07_lp 08_default_posterior; do
  echo "=== $f"; Rscript $f.R 2>&1 | grep -v "tidyselect\|^ℹ\|^  # \|^  data\|^$\|^See <\|^Warning message"
done

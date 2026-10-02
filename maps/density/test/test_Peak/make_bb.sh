#!/usr/bin/env bash
# Rebuilds the bigBed fixtures used by test_Peak.py from the BED files in this directory.
set -euo pipefail
cd "$(dirname "$0")"
for bed in *.bed; do
    bedToBigBed -type=bed6 "$bed" chr1.sizes "$bed.bb"
done

#!/usr/bin/env bash
# Build the alfak2 flowcharts (PDF) included by doc/alfak2_method.tex.
#   doc/flowcharts/build.sh
set -euo pipefail
cd "$(dirname "$0")"
for chart in flow_overview flow_core flow_anchored flow_pairwise flow_state_space; do
  pdflatex -interaction=nonstopmode -halt-on-error "${chart}.tex" > /dev/null
  echo "built ${chart}.pdf"
done
rm -f ./*.aux ./*.log

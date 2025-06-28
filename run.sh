#!/usr/bin/env bash

TIMESTAMP=$(date +%Y%m%d_%H%M%S)

docker run --rm -v $(pwd):/workspace nf-spatial:latest nextflow run nf_spatial.nf \
  -with-report reports/report_${TIMESTAMP}.html \
  -with-trace reports/trace_${TIMESTAMP}.txt \
  -with-dag reports/flowchart_${TIMESTAMP}.png \
  -with-timeline reports/timeline_${TIMESTAMP}.html
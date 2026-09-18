#!/bin/bash

################################################################################
# CellBender Run
################################################################################

METRICS_FILE="$(pwd)/HF_Brain_166_S2/outs/metrics_summary.csv"

EXPECTED_CELLS=$(python3 - <<EOF
import pandas as pd

df = pd.read_csv("${METRICS_FILE}")
print(int(df["Estimated Number of Cells"].iloc[0]))
EOF
)

echo "Estimated cells detected: ${EXPECTED_CELLS}"

docker run -it \
    -v $(pwd)/HF_Brain_166_S2/outs/multi/count:/data \
    us.gcr.io/broad-dsde-methods/cellbender:latest \
    cellbender remove-background \
    --input /data/raw_feature_bc_matrix.h5 \
    --output /data/cellbender_output.h5 \
    --expected-cells ${EXPECTED_CELLS} \
    --low-count-threshold 20

echo "CellBender completed"

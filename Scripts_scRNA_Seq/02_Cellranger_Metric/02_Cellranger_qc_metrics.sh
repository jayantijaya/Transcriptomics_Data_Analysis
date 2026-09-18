#!/bin/bash
###################### CELLRANGER QC METRICS ################################
# INPUT CSV FILE
METRICS_FILE="$(pwd)/HF_Brain_166_S2/outs/per_sample_outs/HF_Brain_166_S2/metrics_summary.csv"

# OUTPUT FILE
OUTPUT_FILE="$(pwd)/cellranger_metrics.txt"

extract_metric () {
    grep "$1" "$METRICS_FILE" | head -1 | \
    python3 -c "
import csv,sys
row=next(csv.reader(sys.stdin))
print(row[-1])
"
}

{
echo "================================================="
echo "CELLRANGER QC SUMMARY"
echo "================================================="

metrics=(
"Cells"
"Median UMI counts per cell"
"Median genes per cell"
"Total genes detected"
"Mean reads per cell"
"Confidently mapped reads in cells"
"Sequencing saturation"
"Valid UMIs"
"Valid barcodes"
"Q30 RNA read"
"Q30 UMI"
"Q30 barcodes"
"Reads mapped to probe set"
"Reads confidently mapped to probe set"
"Reads half-mapped to probe set"
"Reads split-mapped to probe set"
)

for metric in "${metrics[@]}"
do
    echo "$metric: $(extract_metric "$metric")"
done


} > "$OUTPUT_FILE"

echo "Metrics saved to $OUTPUT_FILE"

#!/bin/bash
################################################################################
# Cellranger multi run
################################################################################
/home/hbp/Softwares/cellranger-9.0.1/bin/cellranger multi \
  --id HF_Brain_166_S2 \
  --csv /home/hbp/Softwares/config.csv \
  --localmem 160
  
echo "cellranger completed"

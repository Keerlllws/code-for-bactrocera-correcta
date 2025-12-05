#!/bin/bash

# ============================================================
# Script Name: run_admixture.sh
# Description: Pipeline for VCF conversion, Admixture analysis, 
#              and visualization.
# ============================================================

# Stop the script if any command fails
set -e

echo ">>> Step 1: Converting VCF format to PLINK BED format..."
# Convert VCF to BED/BIM/FAM
plink --vcf ../00.filter/all.LDfilter.vcf \
      --make-bed \
      --out all \
      --allow-extra-chr \
      --keep-allele-order \
      --set-missing-var-ids @:#

echo ">>> Step 2: Running ADMIXTURE analysis (K=2 to K=4)..."
# Loop from K=2 to K=4
# -j2: Use 2 threads
# --cv: Calculate cross-validation error
for K in {2..4}
do
    echo "Processing K=$K..."
    admixture --cv -j2 all.bed $K > admix.${K}.log 2>&1
done

echo ">>> Step 3: Summarizing results and plotting..."
# Create output directory if it doesn't exist (-p)
mkdir -p result

# Copy ancestry fraction files (.Q) to result directory
cp ./*.Q result/

# Run R script for visualization
# Arguments: [Input Dir] [Fam File] [Sample Order File] [Output Prefix]
Rscript ../script/draw_admixture.R \
        result \
        all.fam \
        ../data/sample.pop.order \
        structure

echo ">>> All analysis finished successfully!"
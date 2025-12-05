#!/bin/bash

# ============================================================
# Script Name: run_ld_decay.sh
# Description: Calculate Linkage Disequilibrium (LD) decay 
#              statistics and generate plots for subpopulations.
# Tools Required: PopLDdecay, Plot_OnePop.pl, Plot_MultiPop.pl
# ============================================================

# Stop the script if any command fails
set -e

# ------------------------------------------------------------
# Step 1: Calculate LD (r^2) statistics for subpopulations
# ------------------------------------------------------------
echo ">>> Step 1: Calculating LD statistics for SC and YZR subpopulations..."

# Run for SC population
# -InVCF: Input VCF file
# -SubPop: List of samples in the subpopulation
# -MaxDist: Max distance (kb) between SNP pairs to calculate
# -OutStat: Output prefix
PopLDdecay -InVCF ../data/all.vcf \
           -SubPop ../data/pop.SC.table \
           -MaxDist 500 \
           -OutStat pop.SC.stat

# Run for YZR population
PopLDdecay -InVCF ../data/all.vcf \
           -SubPop ../data/pop.YZR.table \
           -MaxDist 500 \
           -OutStat pop.YZR.stat

# Note: The output files will be automatically gzipped (e.g., pop.SC.stat.gz)

# ------------------------------------------------------------
# Step 2: Plot LD decay for individual populations
# ------------------------------------------------------------
echo ">>> Step 2: Generating individual plots..."

# Plot for SC
Plot_OnePop.pl -inFile pop.SC.stat.gz -output pop.SC.ld

# Plot for YZR
Plot_OnePop.pl -inFile pop.YZR.stat.gz -output pop.YZR.ld

# ------------------------------------------------------------
# Step 3: Plot LD decay for multiple populations together
# ------------------------------------------------------------
echo ">>> Step 3: Generating combined multi-population plot..."

# 3.1 Prepare the configuration file
# This command lists the stat files and extracts the population name (2nd field)
# Output format: [filename] [population_name]
ls pop.*.stat.gz | awk -F"." '{ print $0"\t"$2 }' > ld_stat.list

echo "Configuration file (ld_stat.list) generated:"
cat ld_stat.list

# 3.2 Generate the multi-population plot
# -keepR: Keep the intermediate R script used for plotting
Plot_MultiPop.pl -inList ld_stat.list \
                 -output ld_stat.multi \
                 -keepR

echo ">>> LD Decay analysis finished successfully!"
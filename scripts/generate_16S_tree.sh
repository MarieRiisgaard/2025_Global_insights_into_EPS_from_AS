#!/bin/bash

# Load Conda environment
conda activate tree_env

# MIDAS database
TAXDB="/databases/midas/MiDAS5.3_20240320/FLASVs.fa"
      #head $TAXDB

# Input directory with all ASV files (modify as needed)
ASV_DIR="data/subset_midas_database"

# Output base directory
RESULTS_BASE="results/make_tree_file"

# Check if the output base directory exists; create it if it doesn't
mkdir -p $RESULTS_BASE

# Process each keep_seq file in the ASV_DIR
for keep_seq in $ASV_DIR/*.txt; do
  # Extract the base name of the keep_seq file (e.g., 20240402_MIDAS_52_taxa_subset_top_15_genera_to_tree)
  BASENAME=$(basename "$keep_seq" .txt)

  # Create output directories for the current file
  OUTPUT_DIR="$RESULTS_BASE/$BASENAME"
  mkdir -p $OUTPUT_DIR

  # Define file paths for outputs
  midas_sub="$OUTPUT_DIR/${BASENAME}_MIDAS_subset.fa"
  mafft_msa="$OUTPUT_DIR/${BASENAME}_mafft_output.fasta"
  msa_trimmed="/home/bio.aau.dk/hd95lp/Projects/RETHINK/global_EPS_survay/$OUTPUT_DIR/${BASENAME}_mafft_output_trimmed.msa"
  IQ_tree_out="$OUTPUT_DIR/${BASENAME}_IQ_tree_output.tree"

  echo "Processing file: $keep_seq"
  echo "Results will be stored in: $OUTPUT_DIR"

  # Step 1: Subset MIDAS database
  module load seqtk
  seqtk subseq $TAXDB $keep_seq > $midas_sub
  echo "Created MIDAS subset: $midas_sub"

  # Step 2: Perform MSA using MAFFT
  module load MAFFT/7.402-foss-2018a-with-extensions
  mafft --auto --op 2 $midas_sub > $mafft_msa
  echo "Created MSA: $mafft_msa"

  # Step 3: Trim MSA using TrimAl
  module purge
  module load TrimAl/1.4.1-foss-2020b
  trimal -in $mafft_msa -out $msa_trimmed -gt 0.2
  echo "Trimmed MSA: $msa_trimmed"

#OUTPUT_DIR=results/make_tree_file/20241126_MIDAS53_genera_min_0.1
#msa_trimmed=/home/bio.aau.dk/hd95lp/Projects/RETHINK/global_EPS_survay/results/make_tree_file/20241126_MIDAS53_genera_min_0.1/20241126_MIDAS53_genera_min_0.1_mafft_output_trimmed.msa
#IQ_tree_out="$OUTPUT_DIR/20241126_MIDAS53_genera_min_0.1_IQ_tree_output.tree"

  # Step 4: Build phylogenetic tree using IQ-TREE
  echo "Building tree for: $msa_trimmed"
  cd $OUTPUT_DIR
  iqtree2 -s $msa_trimmed -m SYM+I+R10 -B 1000 -T AUTO --threads-max 20
  echo "Tree saved to: $IQ_tree_out"
done

# Deactivate Conda environment
conda deactivate

echo "Processing completed for all keep_seq files."


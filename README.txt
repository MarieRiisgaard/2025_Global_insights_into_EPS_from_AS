# Code for: *Global insights into extracellular polymeric substances from activated sludge: Yield, composition, and microbial communities*

This repository contains the scripts and data processing workflows used in the analysis presented in the paper:

*Global insights into extracellular polymeric substances from activated sludge: Yield, composition, and microbial communities*  
(*Journal name, year, DOI – add when available*)

---

## 📁 Repository Structure

- `scripts/` – R and shell scripts used for data analysis and figure generation  
- `data/` – Input data used in the study (or instructions to download, if large)  
- `output/` – Plots and files used in the publication  
- `env/` – Conda environment YAML files to reproduce the software environment

---

## 📊 Main Analyses

- Amplicon sequencing data processing using modified [nanopore_16Samp](https://github.com/cmc-aau/nanopore_16Samp) pipeline  
- Raw sequnecing data availible at NCBI bioproject PRJNA1255528 (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1255528)
- Taxonomic classification with MiDAS 5.3 down to genus level  
- Diversity and abundance analysis 
- Phylogenetic tree construction using IQ-TREE  
- Analysis of the extracted EPS from the across global activated sludge samples including: 
	- EPS yield analysis 
	- EPS protein and polysaccharide content
	- FTIR spectra and PCA of the EPS

---

## ⚙️ Environment

To reproduce the analyses, create the Conda environment:

```bash
conda env create -f env/tree_env.yaml
conda activate tree_env

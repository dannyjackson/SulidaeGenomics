# SulidaeGenomics
 Code used to analyze phylogenetics, introgression, and differentiation between booby species

### Contents
This repository contains code to analyze medium coverge whole-genome sequence data from 7 species of boobies (family Sulidae, genus Sula). The steps involved in the analysis are listed below. All scripts in this repository begin with the number associated with the step described below.

1. Download the raw fasta files from the Sequence Read Archive
2. Trim these fastas and assess the quality of the resulting files
3. Align these fastas to the great cormorant (Phalacrocorax carbo) genome (NCBI RefSeq assembly GCF_963921805.1) and clean up the resulting bam files
4. Analyze phylogenetic relationships between and within species
5. Test for introgression between hybridizing pairs of species
6. In instances with a significant test for introgression, scan the genomes for evidence of introgressed regions
7. Use genome scans to identify putative regions of the genome associated with species-specific traits

### Paper status
All codes used in this repository were written to analyze the whole genome sequence data for this paper: https://www.authorea.com/doi/full/10.22541/au.165354126.61061627

Based on feedback from reviewers, we are repeating our analyses using genotype likelihoods instead of genotype calls. This repository is in the process of being updated to reflect that change as we work through the updated pipeline. Additionally, a new chromosome level reference genome was published in the time since we started this project and we are excited about the additional resolution that this provides. This has enabled us to incorporate new sliding window analyses into the next iteration of the paper, with exciting preliminary findings. Stay tuned!
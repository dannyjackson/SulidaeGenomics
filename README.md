# Project Overview
This repository contains scripts and code used to analyze medium-coverage whole-genome sequence data from seven Sulidae species to investigate speciation processes.

### Central Questions
To understand speciation processes in the context of gene flow, we sequenced and analyzed 29 short-read booby genomes and assembled a reference northern gannet genome to (1) test for introgression, (2) characterize genomic patterns of divergence across species, and (3) investigate the gene content and evolution of the neo-sex chromosome. 

### Abbreviations:
   - **RFBO / rfbo**: red-footed booby (*Sula sula*)
   - **BRBO / brbo**: brown booby (*Sula leucogaster*)
   - **MABO / mabo**: masked booby (*Sula dactylatra*)
   - **NABO / nabo**: Nazca booby (*Sula granti*)
   - **BFBO / bfbo**: blue-footed booby (*Sula nebouxii*)
   - **PEBO / pebo**: Peruvian booby (*Sula variegata*)
   - **COBO / cobo**: Cocos booby (*Sula brewsteri*)

### Pipeline
The pipeline investigates 5 major themes, each assigned to a subdirectory of scripts. Within each subdirectory are markdown files, which can be followed in order to replicate the complete pipeline. Each markdown filename follows a structured format:
**[Stage Letter][Step Number]_[Description]**
   - The **letter** represents the pipeline stage.
   - The **number** indicates the order in which the scripts should be excecuted within that stage.
   - The **description** briefly summarizes the purpose of the script.
Each markdown file calls various scripts within this repository, all labeled using a similar syntax.

## Pipeline Stages

0. **Preprocessing**  
1. **Variant Calling**  
2. **Neo-Sex Chromosome Analysis**  
3. **Demographic History Analysis**  
4. **Introgression**  
5. **Divergence**  


All sequence data wille be available on the **Sequence Read Archive** upon publicaton.

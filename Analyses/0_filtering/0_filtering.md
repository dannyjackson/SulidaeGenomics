# Scripts to prepare raw sequence read data for analyses
Steps are as follows:
1. download
2. trim
3. align and sort
4. clean aligned files
5. ?
## Step 1: Download files
If needed, download the fasta files from the Sequence Read Archive

```sbatch 1_download_fastas.sh```

## Step 2: Trim reads
### A. Perform quality control analyses on untrimmed files:
    
```sbatch 2a_pre_fastqc.sh```
    
### B. Trim polyg tails:

```sbatch 2b_polyg_trimming.sh```

### C. Trim on various settings using Trimmomatic:

```sbatch 2c_trimming.sh```

### D. Perform quality control analyses on trimmed files
```sbatch 2d_post_fastqc.sh```


## Step 3: Align and sort
First, download the reference genome from NCBI (can be run interactively):

```
~/programs/datasets download genome accession GCA_031468815.1 --include gff3,rna,cds,protein,genome,seq-report

unzip ncbi_dataset.zip
```
### A. Align fastas to reference genome:

```sbatch --array=1-34%10 3a_align.sh```

### B. Sort bams:

```sbatch --array=1-34%10 3b_sort.sh```

### C. Add read groups and mark duplicates:

```sbatch --array=1-34%10 3c_readgroups_markdups.sh```

Rename sorted bams using species codes and sample numbers, rather than SRR identifiers:
```
conv="/xdisk/mcnew/dannyjackson/sulidae/referencelists/filenameconversion.txt"

while read -r old new; do
  for f in ${old}*; do
    if [[ -e "$f" ]]; then
      newname="${f/$old/$new}"
      mv "$f" "$newname"
    fi
  done
done < "$conv"
```
## Step 4: Clean aligned files

### A. Clip overlap:
```
sbatch --array=1-29%10 4a_clipoverlap.sh
```
### B-D. Indel realignment:

#### First, make a dictionary for the reference genome (necessary for doing indel realignment).
```
REF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna
DICT=/xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.dict

module load picard/2.23.4 

picard CreateSequenceDictionary -R ${REF} -O ${DICT}
```
#### B. Index the overlap-clipped bam files:
```
sbatch --array=1-29%10 4b_indexclip.sh
```
#### C. Make indel maps:
```
sbatch --array=1-29%10 4c_indelmaps.sh
```
#### D. Realign around indels:
```
sbatch --array=1-29%10 4d_indelrealignment.sh
```
Reorganize file locatons:
```
mv /xdisk/mcnew/dannyjackson/sulidae/datafiles/bamfiles/*final* /xdisk/mcnew/dannyjackson/sulidae/datafiles/finalbams/
```

# Step 5: Compute depth of coverage across bams and contigs
A. Compute depth across the genome per individual
```
sbatch --array=1-29%10 5a_bamstats.sh
```
B. Calculate mean depth
```
sbatch --account=mcnew \
--job-name=calcdepthparallel \
--partition=standard \
--mail-type=ALL \
--output=slurm_output/output.calcdepthparallel.%j \
--nodes=1 \
--ntasks-per-node=1 \
--cpus-per-task=12 \
--time=6:00:00 \
calcdepth.parallel.sh
```
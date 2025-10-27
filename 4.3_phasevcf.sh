# Phase data


#!/usr/bin/env bash
#SBATCH --job-name=phasevcf
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --time=24:00:00
#SBATCH --output=slurm_output/phasevcf.%A_%a.out
#SBATCH --mail-type=ALL

# sbatch --array=1-29 phasing.sh /xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplecodes.txt

set -euo pipefail

module load python

LIST=/xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplecodes.txt

IND="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST" | tr -d '\r')"
BAMDIR=/xdisk/mcnew/dannyjackson/sulidae/datafiles/finalbams
BAMFILE=${BAMDIR}/${IND}.final.bam
PROJDIR=/xdisk/mcnew/dannyjackson/sulidae/
OUTDIR=${PROJDIR}/datafiles/snpable_masks
REF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna # reference genome fasta
CHRLIST=${PROJDIR}/referencelists/CONTIGS.txt

echo "working with individual $IND"
for s in `cat ${CHRLIST}`;
    do echo "working with scaffold $s";
    if [ -f ${OUTDIR}/phasedvcfs/${IND}.${s}.whatshap.samtools.vcf.gz ]; then
            echo "phased VCF already exists; moving onto next scaffold";
    else
            echo "phased VCF does not exist; phasing VCF for scaffold $s"
            whatshap phase --reference ${REF} --ignore-read-groups -o ${OUTDIR}/phasedvcfs/${IND}.${s}.whatshap.samtools.vcf.gz ${OUTDIR}/vcf/${IND}.${s}.samtools.vcf ${BAMFILE}
            whatshap stats --tsv=${OUTDIR}/phasedvcfs/stats/${IND}.${s}.${IND}.${s}.minDP10.whatshap.stats.tsv ${OUTDIR}/phasedvcfs/${IND}.${s}.whatshap.samtools.vcf.gz ;
    fi;
done

echo "finished with individual $IND"
date

sbatch --array=1-29 phasing.sh /xdisk/mcnew/dannyjackson/sulidae/referencelists/allsamplecodes.txt

# rename samples in vcfs

echo "working with individual $IND"
for s in `cat ${CHRLIST}`;
    do echo "working with scaffold $s";
            echo "phased VCF does not exist; phasing VCF for scaffold $s"
            whatshap phase --reference ${REF} --ignore-read-groups -o ${OUTDIR}/phasedvcfs/${IND}.${s}.whatshap.samtools.vcf.gz ${OUTDIR}/vcf/${IND}.${s}.samtools.vcf ${BAMFILE}
            whatshap stats --tsv=${OUTDIR}/phasedvcfs/stats/${IND}.${s}.${IND}.${s}.minDP10.whatshap.stats.tsv ${OUTDIR}/phasedvcfs/${IND}.${s}.whatshap.samtools.vcf.gz ;
    fi;
done


#!/usr/bin/env bash
set -euo pipefail

shopt -s nullglob
for vcf in *.whatshap.samtools.vcf.gz; do
  # specimen name is the bit before the first dot in the filename
  base="$(basename "$vcf")"
  newsample="${base%%.*}"
  echo "Renaming sample in $vcf to $newsample"
  # sanity check: expect exactly 1 sample in the VCF
  ns=$(bcftools query -l "$vcf" | wc -l)
  if [ "$ns" -ne 1 ]; then
    echo "Skipping $vcf: has $ns samples (expected 1)" >&2
    continue
  fi

  # replace the single sample name with the new one
  out=renamedvcfs/"$vcf"
  printf "%s\n" "$newsample" > "${out%.vcf.gz}.names"
  bcftools reheader -s "${out%.vcf.gz}.names" -o "$out" "$vcf"
  tabix -f "$out"

  rm -f "${out%.vcf.gz}.names"
  echo "Renamed sample in $vcf -> $newsample"
done

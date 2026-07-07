# Genome annotation with LiftOff/LiftOn

## Define input files

```bash
# Target genome (Brown Booby, haplotype 2)
REF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna

# Source annotation (Great Cormorant)
PhaCarFasta=/xdisk/mcnew/dannyjackson/sulidae/datafiles/other_ref_genomes/great_cormorant/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna
PhaCarGFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/other_ref_genomes/great_cormorant/ncbi_dataset/data/GCF_963921805.1/genomic.gff
```

## Download Great Cormorant reference files

```bash
~/programs/datasets download genome accession GCF_963921805.1 \
    --include genome,gff3,rna,cds,protein,seq-report
```

## Install LiftOff

```bash
conda install -c bioconda liftoff
```

---

# LiftOff annotation transfer

Transfers Great Cormorant gene annotations onto the Brown Booby reference genome.

```bash
#!/usr/bin/env bash
#SBATCH --job-name=LiftOff_PhaCar
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=24
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --output=slurm_output/LiftOff_PhaCar.%A_%a.out
#SBATCH --mail-type=ALL

# Target genome
REF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna

# Source genome and annotations
PhaCarFasta=/xdisk/mcnew/dannyjackson/sulidae/datafiles/other_ref_genomes/great_cormorant/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna
PhaCarGFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/other_ref_genomes/great_cormorant/ncbi_dataset/data/GCF_963921805.1/genomic.gff

source ~/.bashrc
module load micromamba
eval "$(micromamba shell hook --shell bash)"

# Create environment once:
# micromamba create -n LiftGFF -c conda-forge -c bioconda liftoff

micromamba activate LiftGFF

liftoff \
    -g $PhaCarGFF \
    -o /xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/GCA_031468815.1_bMorBas2.PhaCar.hap2_genomic_lifted.gff \
    -p 24 \
    $REF \
    $PhaCarFasta
```

---

# Generate background gene lists

Extract gene IDs from the lifted GFF for downstream enrichment analyses.

### All genes

```bash
grep 'ID=gene' ${GENEFILE} \
    | awk '{OFS="\t"} {split($9,arr,";"); print arr[1]}' \
    | sed 's/ID=gene-//g' \
    | sort -u \
    > ${GENENAMES}
```

### Autosomal genes only

Excludes chromosomes `CM062600.1` and `CM062610.1`.

```bash
awk -F'\t' 'BEGIN{OFS="\t"}
  $1 ~ /^CM/ &&
  $1 !~ /^CM062600\.1$/ &&
  $1 !~ /^CM062610\.1$/ &&
  $9 ~ /ID=gene/ {
    split($9,a,";");
    sub(/^ID=gene-/,"",a[1]);
    print a[1]
}' GCA_031468815.1_bMorBas2.PhaCar.hap2_genomic_lifted.gff \
| sort -u \
> GCA_031468815.1_bMorBas2.PhaCar.hap2_genomic_lifted.backgroundgenes.autosomes.txt
```

---

# LiftOn annotation transfer

Runs LiftOn using the Great Cormorant annotation and reports duplicated gene copies.

```bash
#!/usr/bin/env bash
#SBATCH --job-name=LiftOn_PhaCar
#SBATCH --partition=standard
#SBATCH --account=mcnew
#SBATCH --cpus-per-task=24
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --output=slurm_output/LiftOn_PhaCar.%A_%a.out
#SBATCH --mail-type=ALL

source ~/.bashrc
module load micromamba
eval "$(micromamba shell hook --shell bash)"

# Install once:
# pip install lifton
# micromamba install miniprot

micromamba activate LiftGFF

# Target genome
REF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/reference_genome/ncbi_dataset/data/GCA_031468815.1/GCA_031468815.1_bMorBas2.hap2_genomic.fna

# Source genome and annotations
PhaCarFasta=/xdisk/mcnew/dannyjackson/sulidae/datafiles/other_ref_genomes/great_cormorant/ncbi_dataset/data/GCF_963921805.1/GCF_963921805.1_bPhaCar2.1_genomic.fna
PhaCarGFF=/xdisk/mcnew/dannyjackson/sulidae/datafiles/other_ref_genomes/great_cormorant/ncbi_dataset/data/GCF_963921805.1/genomic.gff

lifton \
    -g $PhaCarGFF \
    -o /xdisk/mcnew/dannyjackson/sulidae/datafiles/liftoff_annotations/GCA_031468815.1_bMorBas2.PhaCar.hap2_genomic_lifted.on.copies.gff \
    -copies \
    -t 24 \
    $REF \
    $PhaCarFasta
```
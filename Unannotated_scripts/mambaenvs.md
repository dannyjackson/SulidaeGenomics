# Creating micromamba environments
# micromamba clean --index-cache -y

###########################
# elgato 
###########################
name: r_elgato
channels:
  - conda-forge
dependencies:
  - r-base=4.4.*
  - r-ggplot2=3.5.*
  - r-tidyr=1.3.*
micromamba create -f /xdisk/mcnew/dannyjackson/sulidae/yamls/r_elgato.yml -y

micromamba activate r_elgato
micromamba install -c bioconda mcscanx -y
micromamba install -c conda-forge r-devtools r-usethis r-miniui r-pkgdown -y
micromamba install -c conda-forge -c bioconda bioconductor-biostrings
micromamba install -c conda-forge -c bioconda bioconductor-GenomicRanges
micromamba install -c conda-forge -c bioconda bioconductor-Rsamtools
micromamba install -c conda-forge r-data.table
micromamba install -c conda-forge r-r.utils
micromamba install -c conda-forge r-readr
micromamba install -c conda-forge r-ggplot2
micromamba install -c conda-forge r-patchwork
micromamba install -c conda-forge r-gggenes
micromamba install -c conda-forge r-ggnewscale
micromamba install -c conda-forge r-cowplot

micromamba create -n genomics_elgato -c bioconda -c conda-forge seqkit

##########################################################################################

name: bcftools_env
channels:
  - bioconda
  - conda-forge
  - defaults
dependencies:
  - bcftools
  - htslib
  - samtools          # optional, for BAM/CRAM work
  - tabix                  # for indexing bgzip files
  - parallel               # handy for batch VCF operations

micromamba create -f bcftools_elgato.yml

micromamba activate bcftools_elgato

export BCFTOOLS_PLUGINS="$CONDA_PREFIX/libexec/bcftools"

########################################################################################
name: msmc_im_elgato
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.11
  - numpy
  - scipy
  - matplotlib
  - pip

micromamba create -f msmc_im_elgato.yml

micromamba install -c conda-forge pandas

micromamba activate msmc_im_elgato

########################################################################################


###########################
# puma 
###########################
name: r_puma
channels:
  - conda-forge
dependencies:
  - r-base=4.4.*
  - r-ggplot2=3.5.*
  - r-tidyr=1.3.*
micromamba create -f /xdisk/mcnew/dannyjackson/sulidae/yamls/r_puma.yml -y

micromamba activate r_puma
micromamba install -c bioconda mcscanx -y
micromamba install -c conda-forge r-devtools r-usethis r-miniui r-pkgdown -y
micromamba install -c conda-forge -c bioconda bioconductor-biostrings
micromamba install -c conda-forge -c bioconda bioconductor-GenomicRanges
micromamba install -c conda-forge -c bioconda bioconductor-Rsamtools
micromamba install -c conda-forge r-data.table
micromamba install -c conda-forge r-dbscan
micromamba install -c conda-forge r-r.utils


###########################
# ocelote 
###########################
name: r_ocelote
channels:
  - conda-forge
dependencies:
  - r-base=4.4.*
  - r-ggplot2=3.5.*
  - r-tidyr=1.3.*

micromamba create -f /xdisk/mcnew/dannyjackson/sulidae/yamls/r_ocelote.yml -y

micromamba activate r_ocelote
micromamba install -c conda-forge r-devtools r-usethis r-miniui r-pkgdown -y
micromamba install -c conda-forge -c bioconda bioconductor-biostrings
micromamba install -c conda-forge -c bioconda bioconductor-GenomicRanges
micromamba install -c conda-forge -c bioconda bioconductor-Rsamtools
micromamba install -c conda-forge r-data.table
micromamba install -c conda-forge r-r.utils
micromamba install -c conda-forge r-readr
micromamba install -c conda-forge r-dbscan
micromamba install -c conda-forge r-igraph
micromamba install -c conda-forge r-tidyverse
micromamba install -c conda-forge r-qqman r-hexbin r-ggrepel r-dplyr r-RColorBrewer
micromamba install -c conda-forge r-patchwork
micromamba install -c conda-forge r-genespace
micromamba install orthofinder
micromamba install -c conda-forge biopython
micromamba install -c conda-forge -c bioconda bioconductor-Gviz
micromamba install -c conda-forge -c bioconda bioconductor-rtracklayer
micromamba install -c conda-forge r-ape
micromamba install -c conda-forge r-phangorn
micromamba install -c bioconda r-ggtree
micromamba install r-phytools

micromamba install -c conda-forge numpy
micromamba install -c conda-forge pandas
micromamba install -c conda-forge matplotlib


###########################
# ocelote  -- pcangsd
###########################
# install pcangsd

micromamba create -n pcangsd \
  -c conda-forge -c bioconda \
  pcangsd

micromamba activate pcangsd


micromamba install -c conda-forge r-base=4.4.*
micromamba install -c conda-forge r-ggplot2=3.5.*
micromamba install -c conda-forge r-tidyr=1.3.*
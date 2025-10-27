# Installing ggtree


module load micromamba

micromamba create -n r42_ggtree -c conda-forge -c bioconda \
  r-base=4.2 \
  r-ggplot2 \
  r-ggfun \
  bioconductor-ggtree \
  bioconductor-treeio \
  bioconductor-tidytree \
  -y

module purge
# make sure your CentOS7 R-4.2 user lib is first in .libPaths():
# echo 'R_LIBS_USER=~/R/library_4.2_cent7' > ~/.Renviron
# mkdir -p ~/R/library_4.2_cent7
# module load R/4.2
# echo 'R_LIBS_USER=~/R/library_4.4_cent7' > ~/.Renviron
echo 'R_LIBS_USER=~/R/library_4.4_rocky9' > ~/.Renviron
mkdir -p ~/R/library_4.4_rocky9
module load R/4.4.0

# Keep everything source and avoid “helpful” binary pulls
Sys.unsetenv("RSPM")
options(repos = c(CRAN = "https://cloud.r-project.org"), pkgType = "source")
Sys.setenv(MAKEFLAGS = "-j8")

# 1) Bring ggplot2 and its light deps up to date (no system libs needed)
install.packages(c(
  "rlang","vctrs","glue","digest","withr","isoband","farver","scales","gtable",
  "pillar","lifecycle"  # safe, pure R/C++ deps
))
install.packages("ggplot2")   # >= 3.4.x provides check_linewidth()

# 2) Bioc side: install ggtree deps without updating curl/openssl/etc.
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(version = "3.16", ask = FALSE, update = FALSE)

BiocManager::install(c("treeio","tidytree"), ask = FALSE, update = FALSE)
install.packages("ggfun")  # CRAN helper used by ggtree

# 3) Finally, ggtree (do NOT set update=TRUE)
BiocManager::install("ggtree", ask = FALSE, update = FALSE)

# didn't work

remotes::install_github("YuLab-SMU/ggtree")
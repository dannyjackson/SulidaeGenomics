# ld_decay

cd /xdisk/mcnew/dannyjackson/sulidae/analyses/linkagedecay

VCF=/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/chroms/CM062567.1.qualitysort_filtered_mind2.vcf

plink --vcf $VCF --double-id --allow-extra-chr \
--set-missing-var-ids @:# \
--thin 0.1 -r2 gz --ld-window 100 --ld-window-kb 1000 \
--ld-window-r2 0 \
--make-bed --out chr1

# Calculating average LD across set distances

scripts/ld_decay_calc.py -i chr1.ld.gz -o chr1




library(tidyverse)

# set path
my_bins <- "chr1.ld_decay_bins"

# read in data
ld_bins <- read_tsv(my_bins)

# plot LD decay
ggplot(ld_bins, aes(distance, avg_R2)) + geom_line() +
  xlab("Distance (bp)") + ylab(expression(italic(r)^2))

ggsave("ld_decay_plot.png", width = 7, height = 5, dpi = 200)

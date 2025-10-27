# Identify sex of individuals based on bamfiles post-filtering
PROJDIR=/xdisk/mcnew/dannyjackson/sulidae/
OUTDIR=${PROJDIR}/datafiles/bamstats/sexing_by_depth_filtered
BAMDIR=/xdisk/mcnew/dannyjackson/sulidae/datafiles/finalbams/


module load micromamba 

micromamba activate r_puma
mkdir -p ${OUTDIR}
while read -r bird; do
samtools idxstats ${BAMDIR}/"$bird".final.bam | grep 'CM' > ${OUTDIR}/"$bird".idxstats
done <  ${PROJDIR}/referencelists/allsamplecodes.txt 


ls "$PWD"/*idxstats > ${PROJDIR}/referencelists/idxindex.txt


R

library(tidyverse)
library(data.table)
library(ggplot2)
library(readr)
library(dplyr)
library(purrr)
library(stringr)

# list of idxstats files
file_list <- readLines("/xdisk/mcnew/dannyjackson/sulidae/referencelists/idxindex.txt")
sex_file  <- "/xdisk/mcnew/dannyjackson/sulidae/referencelists/sexID.csv"
conv_file <- "/xdisk/mcnew/dannyjackson/sulidae/referencelists/chromconversion.txt"

# read and combine into one dataframe
all_data <- file_list %>%
  set_names(basename(.)) %>%   # name each dataset by filename
  map_dfr(~ read.table(.x, header = FALSE, sep = "\t", stringsAsFactors = FALSE),
          .id = "individual")

sex <- read_csv(sex_file, col_names = c("individual", "sex"))
conv <- read_csv(conv_file, col_names = c("chrom_num", "CHROM"))

# add column names
colnames(all_data) <- c("individual", "chrom", "length", "mapped", "unmapped")

# NORMALIZE DEPTH
# compute per-chrom depth
all_data <- all_data %>%
  mutate(depth = mapped / length)

# calculate normalization factors for each individual
norm_factors <- all_data %>%
  filter(chrom %in% c("CM062567.1", "CM062568.1", "CM062569.1")) %>%
  group_by(individual) %>%
  summarize(norm_factor = mean(depth, na.rm = TRUE))

# join and normalize
all_data <- all_data %>%
  left_join(norm_factors, by = "individual") %>%
  mutate(norm_depth = depth / norm_factor)

all_data <- all_data %>%
  mutate(individual = str_remove(individual, "\\.idxstats$"))

# --- Join with sex and chr ---
depth_sex <- all_data %>%
  left_join(sex, by = "individual") %>%                    # add sex
  left_join(conv, by = c("chrom" = "CHROM")) %>%           # add chromosome number
  mutate(
    sex = str_squish(sex),                                 # trim any stray spaces
    sex = coalesce(sex, "Unknown")                         # replace NA with "Unknown"
  )

# for ordering chromosomes
depth_sex <- depth_sex %>%
  mutate(chrom = factor(chrom, levels = unique(chrom)))

# keep_chroms <- c("CM062567.1", "CM062568.1", "CM062569.1", "CM062600.1")

#depth_sex <- depth_sex %>%
#  filter(chrom %in% keep_chroms) %>%
#  mutate(chrom = factor(chrom, levels = keep_chroms))

drop_chroms <- c("CM062610.1")
depth_sex <- depth_sex %>%
  filter(!chrom %in% drop_chroms) %>%
  mutate(chrom = factor(chrom, levels = unique(chrom)))


# plot
ggplot(depth_sex, aes(x = chrom_num, y = norm_depth, color = individual)) +
  geom_point(position = position_jitter(width = 0.2, height = 0)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(x = "Chromosome", y = "Normalized depth")

ggsave("normalized_depth_plot_ind.png", width = 12, height = 6)



# Ensure chrom is character
depth_sex <- depth_sex %>% mutate(chrom_num = as.character(chrom_num))


# Compute desired factor levels once
lev_base   <- str_sort(unique(depth_sex$chrom_num), numeric = TRUE)
specials   <- c("Z","W","X","MT","M","Un","Unplaced")
lev_order  <- c(setdiff(lev_base, specials),
                specials[specials %in% lev_base])   # push specials to the end (if present)


# Apply levels
depth_sex <- depth_sex %>%
  mutate(
    chrom_num = factor(chrom_num, levels = lev_order),
    sex = factor(sex, levels = c("Female", "Male", "Unknown"))
  )
  
# --- Plot: chromosomes on x-axis, violins by sex ---

pd <- position_dodge(width = 0.8)

gg <- ggplot(depth_sex, aes(x = chrom_num, y = norm_depth, fill = sex)) +
  geom_violin(position = pd, width = 0.75, trim = FALSE, scale = "width", color = NA) +
  stat_summary(fun = median, geom = "point", position = pd, shape = 95, size = 5, color = "black") +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
             alpha = 0.25, size = 0.6, stroke = 0) +
  scale_fill_manual(values = c(Female="#E26D5A", Male="#4F7CAC", Unknown="grey70")) +
  labs(
    x = "Chromosome",
    y = expression(paste("Normalized Depth"))
  ) +
  theme_bw(base_size = 18) +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12)
  )

ggsave("normalized_depth_plot.sex.pdf", width = 12, height = 3)

depth_sex <- depth_sex %>%
  mutate(species = str_extract(individual, "^[A-Z]{4}"))


keep_chroms <- c("1", "2", "3", "29", "33", "Z")

depth_species <- depth_sex %>%
  filter(chrom_num %in% keep_chroms) %>%
  mutate(chrom_num = factor(chrom_num, levels = keep_chroms))

gg <- ggplot(depth_species, aes(x = chrom_num, y = norm_depth, fill = species)) +
  geom_violin(position = pd, width = 0.75, trim = FALSE, scale = "width", color = NA) +
  stat_summary(fun = median, geom = "point", position = pd, shape = 95, size = 5, color = "black") +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
             alpha = 0.25, size = 0.6, stroke = 0) +
  scale_fill_manual(values = c(BFBO="#4FE98C", PEBO="#2894FF", RFBO="#F3447D", BRBO="#000000", MABO="#E8CA3D", NABO="#B32904")) +
  labs(
    x = "Chromosome",
    y = expression(paste("Normalized Depth"))
  ) +
  theme_bw(base_size = 18) +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )

ggsave("normalized_depth_plot.species.pdf", width = 6, height = 3)

# MABO305 is the outlier on chr 33 with high depth

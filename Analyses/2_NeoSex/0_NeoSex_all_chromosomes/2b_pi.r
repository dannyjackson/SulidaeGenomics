# Load packages
library(tidyverse)

# --- Input files ---
pi_file   <- "mean_pi_matrix.tsv"
sex_file  <- "/xdisk/mcnew/dannyjackson/sulidae/referencelists/sexID.csv"
conv_file <- "/xdisk/mcnew/dannyjackson/sulidae/referencelists/chromconversion.txt"

# --- Read data ---
pi  <- read_tsv(pi_file, na = c("", "NA", ".", "NaN"))
sex <- read_csv(sex_file, col_names = c("individual", "sex"))
conv <- read_csv(conv_file, col_names = c("chrom_num", "CHROM"))

# --- Reshape & join ---
pi_sex <- pi %>%
  pivot_longer(-CHROM, names_to = "individual", values_to = "pi", values_drop_na = TRUE) %>%
  left_join(sex, by = "individual") %>%
  left_join(conv, by = "CHROM") %>%       # attach chromosome number
  mutate(sex = replace_na(sex, "Unknown")) %>%
  filter(is.finite(pi))

# Ensure chrom is character
pi_sex <- pi_sex %>% mutate(chrom_num = as.character(chrom_num))


# Compute desired factor levels once
lev_base   <- str_sort(unique(pi_sex$chrom_num), numeric = TRUE)
specials   <- c("Z","W","X","MT","M","Un","Unplaced")
lev_order  <- c(setdiff(lev_base, specials),
                specials[specials %in% lev_base])   # push specials to the end (if present)


# Apply levels
pi_sex <- pi_sex %>%
  mutate(
    chrom_num = factor(chrom_num, levels = lev_order),
    sex = factor(sex, levels = c("Female", "Male", "Unknown"))
  )
  

# --- Plot: chromosomes on x-axis, violins by sex ---
pd <- position_dodge(width = 0.8)

gg <- ggplot(pi_sex, aes(x = chrom_num, y = pi, fill = sex)) +
  geom_violin(position = pd, width = 0.75, trim = FALSE, scale = "width", color = NA) +
  stat_summary(fun = median, geom = "point", position = pd, shape = 95, size = 5, color = "black") +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
             alpha = 0.25, size = 0.6, stroke = 0) +
  scale_fill_manual(values = c(Female="#E26D5A", Male="#4F7CAC", Unknown="grey70")) +
  labs(
    x = "Chromosome",
    y = expression(paste("Heterozygosity (", pi, ")"))
  ) +
  theme_bw(base_size = 18) +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12)
  )


# --- Save ---
ggsave("heterozygosity_by_chromosome.sex.pdf", gg, width = 12, height = 3, dpi = 300)


pi_sex <- pi_sex %>%
  mutate(species = str_extract(individual, "^[A-Z]{4}"))


keep_chroms <- c("1", "2", "3", "29", "33", "Z")

pi_species <- pi_sex %>%
  filter(chrom_num %in% keep_chroms) %>%
  mutate(chrom_num = factor(chrom_num, levels = keep_chroms))

gg <- ggplot(pi_species, aes(x = chrom_num, y = pi, fill = species)) +
  geom_violin(position = pd, width = 0.75, trim = FALSE, scale = "width", color = NA) +
  stat_summary(fun = median, geom = "point", position = pd, shape = 95, size = 5, color = "black") +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
             alpha = 0.25, size = 0.6, stroke = 0) +
  scale_fill_manual(values = c(BFBO="#4FE98C", PEBO="#2894FF", RFBO="#F3447D", BRBO="#000000", MABO="#E8CA3D", NABO="#B32904")) +
  labs(
    x = "Chromosome",
    y = expression(paste("Heterozygosity (", pi, ")"))
  ) +
  theme_bw(base_size = 18) +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12)
  )


# --- Save ---
ggsave("heterozygosity_by_chromosome.species.pdf", gg, width = 6, height = 3, dpi = 300)


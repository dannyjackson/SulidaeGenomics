# Load packages
library(tidyverse)

# --- Input files ---
het_file   <- "all_het_matrix.tsv"
sex_file  <- "/xdisk/mcnew/dannyjackson/sulidae/referencelists/sexID.csv"
conv_file <- "/xdisk/mcnew/dannyjackson/sulidae/referencelists/chromconversion.txt"


# --- Read data ---
het  <- read_tsv(het_file, na = c("", "NA", ".", "NaN"), show_col_types = FALSE)
sex  <- read_csv(sex_file, col_names = c("individual", "sex"), show_col_types = FALSE)
conv <- read_csv(conv_file, col_names = c("chrom_num", "CHROM"), show_col_types = FALSE)

het <- het %>%
  select(CHROM, IND, H_O)

# --- Clean/standardize keys ---
sex <- sex %>%
  mutate(
    individual = str_trim(individual),
    sex = case_when(
      str_to_lower(sex) %in% c("f","female") ~ "Female",
      str_to_lower(sex) %in% c("m","male")   ~ "Male",
      TRUE ~ "Unknown"
    )
  )

conv <- conv %>%
  mutate(CHROM = str_trim(CHROM)) %>%
  distinct(CHROM, .keep_all = TRUE)


het_full <- het %>%
  left_join(sex,  by = c("IND" = "individual")) %>%  # add sex column
  left_join(conv, by = "CHROM")                      # add chrom_num



# Compute desired factor levels once
lev_base   <- str_sort(unique(het_full$chrom_num), numeric = TRUE)
specials   <- c("Z","W","X","MT","M","Un","Unplaced")
lev_order  <- c(setdiff(lev_base, specials),
                specials[specials %in% lev_base])   # push specials to the end (if present)


# Apply levels
het_full <- het_full %>%
  mutate(
    chrom_num = factor(chrom_num, levels = lev_order),
    sex = factor(sex, levels = c("Female", "Male", "Unknown"))
  )
  

# --- Plot: chromosomes on x-axis, violins by sex ---
pd <- position_dodge(width = 0.8)

gg <- ggplot(het_full, aes(x = chrom_num, y = H_O, fill = sex)) +
  geom_violin(position = pd, width = 0.75, trim = FALSE, scale = "width", color = NA) +
  stat_summary(fun = median, geom = "point", position = pd, shape = 95, size = 5, color = "black") +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
             alpha = 0.25, size = 0.6, stroke = 0) +
  scale_fill_manual(values = c(Female="#E26D5A", Male="#4F7CAC", Unknown="grey70")) +
  labs(
    x = "Chromosome",
    y = NULL
  ) +
  theme_bw(base_size = 18) +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
    legend.position = "none"
  )


# --- Save ---
ggsave("heterozygosity_by_chromosome.sex.pdf", gg, width = 12, height = 3, dpi = 300)


het_full <- het_full %>%
  mutate(species = str_extract(IND, "^[A-Z]{4}"))


keep_chroms <- c("1", "2", "3", "6", "17", "24", "29", "33", "Z")

het_species <- het_full %>%
  filter(chrom_num %in% keep_chroms) %>%
  mutate(chrom_num = factor(chrom_num, levels = keep_chroms))

gg <- ggplot(het_species, aes(x = chrom_num, y = H_O, fill = species)) +
  geom_violin(position = pd, width = 0.75, trim = FALSE, scale = "width", color = NA) +
  stat_summary(fun = median, geom = "point", position = pd, shape = 95, size = 5, color = "black") +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
             alpha = 0.25, size = 0.6, stroke = 0) +
  scale_fill_manual(values = c(BFBO="#4FE98C", PEBO="#2894FF", RFBO="#F3447D", BRBO="#000000", MABO="#E8CA3D", NABO="#B32904")) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 18) +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
    legend.position = "none"
  )


# --- Save ---
ggsave("heterozygosity_by_chromosome.species.pdf", gg, width = 6, height = 3, dpi = 300)

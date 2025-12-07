# --- packages ---
library(readr); library(dplyr); library(tidyr); library(stringr); library(ggplot2)

WIN <- 1000000
OUT <- "hets_CM062600.1"   # must match the bash OUT above (or set dynamically)

samples <- readLines(paste0(OUT, ".samples.txt"))
colnames_geno <- c("CHROM","POS", samples)

# Read genotype table: POS numeric; GT columns character
geno <- read_tsv(paste0(OUT, ".genotypes.tsv"),
                 col_names = colnames_geno,
                 col_types = cols(
                   CHROM = col_character(),
                   POS   = col_double(),
                   .default = col_character()
                 ),
                 progress = FALSE)

# Long format
long <- geno |>
  pivot_longer(cols = all_of(samples),
               names_to = "individual",
               values_to = "GT")

# Heterozygote flag (accepts / or |; treats ./. as NA)
is_het <- function(gt) {
  ifelse(is.na(gt) | gt %in% c("./.", ".|."), NA,
         {
           m <- str_match(gt, "^([0-9.]+)[/|]([0-9.]+)$")
           ifelse(is.na(m[,1]), NA,
                  (m[,2] != m[,3]) & (m[,2] != ".") & (m[,3] != "."))
         })
}

# Bin positions to 1 Mbp windows
hets <- long |>
  mutate(
    het = is_het(GT),
    win_start = floor((POS - 1) / WIN) * WIN,
    win_mid   = win_start + WIN/2
  ) |>
  group_by(individual, CHROM, win_start, win_mid) |>
  summarise(
    n_sites = sum(!is.na(het)),
    n_het   = sum(het, na.rm = TRUE),
    het_rate = n_het / n_sites,
    .groups = "drop"
  )

# --- Species assignment ---
# If IDs start with species code (e.g., "BFBO501"), use substring:
hets <- hets |>
  mutate(species = substr(individual, 1, 4))
# If you have a metadata table instead, do something like:
# meta <- read_tsv("samples_metadata.tsv")  # columns: individual, species
# hets <- hets |> left_join(meta, by="individual")

# Optional: filter out low-data windows (helps noisy tails)
min_sites <- 100                         # tweak to your data density
hets_f <- hets |>
  filter(n_sites >= min_sites)

# --- Build a per-species palette mapping: individual -> color ---
# Use a strong qualitative palette; repeat per species (that's fine)
pal_for <- function(n) grDevices::hcl.colors(n, palette = "Dynamic")

inds_by_sp <- hets_f |>
  distinct(species, individual) |>
  arrange(species, individual) |>
  group_by(species) |>
  summarise(ind_list = list(individual), .groups = "drop")

# one palette *per* species; same colors can recur in a different species facet
pal_tbl <- inds_by_sp |>
  rowwise() |>
  mutate(cols = list(setNames(pal_for(length(ind_list)), ind_list))) |>
  ungroup()

pal_individual <- unlist(pal_tbl$cols)  # named vector: names = individual, values = hex colors

sex_file  <- "/xdisk/mcnew/dannyjackson/sulidae/referencelists/sexID.csv"

sex <- read_csv(sex_file, col_names = c("individual", "sex"))

# --- Join with sex and chr ---
hets_sex <- hets_f %>%
  left_join(sex, by = "individual") %>%                    # add sex
  mutate(
    sex = str_squish(sex),                                 # trim any stray spaces
    sex = coalesce(sex, "Unknown")                         # replace NA with "Unknown"
  )


# Apply levels
hets_sex <- hets_sex %>%
  mutate(
    sex = factor(sex, levels = c("Female", "Male", "Unknown"))
  )

# --- Plot ---
p <- ggplot(hets_sex,
            aes(x = win_mid, y = het_rate, color = sex, group = individual)) +
  geom_line(linewidth = 0.2, alpha = 0.95) +
  scale_color_manual(values = c(Female="#E26D5A", Male="#4F7CAC", Unknown="grey70")) +
  labs(
    title = paste0("Z Chromosome (1 Mb bins)"),
    x = "Position (bp)", y = "Heterozygosity", color = "Individual"
  ) +
  theme_bw(base_size = 18) +
  theme(
    panel.spacing.y = unit(0.5, "lines"),
    legend.position = "right",
    legend.key.height = unit(0.35, "lines"),
    legend.text = element_text(size = 18)
  )


ggsave(paste0(OUT, ".hets_1Mb.pdf"), p, width = 12, height = 3, dpi = 300)

library(patchwork)

plots_facets <- hets_f %>%
  split(.$species) %>%
  lapply(function(df) {
    inds <- sort(unique(df$individual))
    pal  <- setNames(pal_for(length(inds)), inds)

    ggplot(df, aes(x = win_mid, y = het_rate,
                   color = individual, group = individual)) +
      geom_line(linewidth = 0.6, alpha = 0.9) +
      scale_color_manual(values = pal) +
      labs(title = unique(df$species), x = "Position (bp)", y = NULL, color = NULL) +
      theme_bw(base_size = 10) +
      theme(legend.position = "right",
            legend.key.height = unit(0.4, "lines"),
            plot.margin = margin(2, 5, 2, 5))
  })

wrap_plots(plots_facets, ncol = 1)
ggsave(paste0(OUT, ".hets_1Mb.faceted.pdf"), width = 12, height = 8, dpi = 300)

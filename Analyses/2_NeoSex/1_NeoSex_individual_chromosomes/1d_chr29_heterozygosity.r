# --- packages ---
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(gggenes)   # for gene arrows or blocks
library(ggnewscale)
library(cowplot)

WIN <- 50000
OUT <- "hets_CM062595.1"   # must match the bash OUT above (or set dynamically)

samples <- sub("\\.bam$", "", basename(readLines(paste0(OUT, ".samples.txt"))))
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

# Bin positions to 50 kb windows
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


# Read GFF PC
genesPC <- read_tsv("/xdisk/mcnew/dannyjackson/sulidae/analyses/neosex_inference/genes_chr29/chr29.PC.genes.gff", comment = "#", col_names = FALSE)

colnames(genesPC) <- c("seqid", "source", "type", "start", "end", "score",
                     "strand", "phase", "attributes")

# Filter only genes
genesPC <- genesPC %>% filter(type == "gene")

# Extract gene name from attributes
genesPC <- genesPC %>%
  mutate(gene_name = sub(".*Name=([^;]+);.*", "\\1", attributes))

# assign stacking level
# Sort and assign "y position" to stack genes
genesPC <- genesPC %>%
  arrange(start) %>%
  mutate(
    y = case_when(
      strand == "+" ~ 1,
      strand == "-" ~ 2,
      TRUE ~ NA_real_     # optional: handles missing strand values
    )
  )


# Read GFF MB

genesMB <- read_tsv("/xdisk/mcnew/dannyjackson/sulidae/analyses/neosex_inference/genes_chr29/chr29.MB.genes.gff", comment = "#", col_names = FALSE)

colnames(genesMB) <- c("seqid", "source", "type", "start", "end", "score",
                     "strand", "phase", "attributes")

# Filter only genes
genesMB <- genesMB %>% filter(type == "gene")

# Extract gene name from attributes
genesMB <- genesMB %>%
  mutate(gene_name = sub(".*Name=([^;]+);.*", "\\1", attributes))

# assign stacking level
# Sort and assign "y position" to stack genes
genesMB <- genesMB %>%
  arrange(start) %>%
  mutate(
    y = case_when(
      strand == "+" ~ 1,
      strand == "-" ~ 2,
      TRUE ~ NA_real_     # optional: handles missing strand values
    )
  )


# Read GFF RE
genesRE <- read_tsv("/xdisk/mcnew/dannyjackson/sulidae/analyses/neosex_inference/genes_chr29/chr29.RE.gff", comment = "#", col_names = FALSE)

colnames(genesRE) <- c("seqid", "source", "type", "start", "end", "score",
                     "strand", "phase", "attributes")

# Filter only genes
# genesRE <- genesRE %>% filter(type == "gene")

# Extract gene name from attributes
genesRE <- genesRE %>%
  mutate(gene_name = sub(".*Name=([^;]+);.*", "\\1", attributes))



# combine with heterozygosity plot:


xmax <- 4730139


# Heterozygosity plot (your code)
p1 <- ggplot(hets_sex,
             aes(x = win_mid, y = het_rate, color = sex, group = individual)) +
  geom_line(linewidth = 0.2, alpha = 0.95) +
  scale_color_manual(values = c(Female="#E26D5A", Male="#4F7CAC", Unknown="grey70")) +
  labs(
    x = NULL, y ="H",
  ) +
  scale_x_continuous(limits = c(0, xmax), expand = c(0, 0)) +
  theme_bw(base_size = 15) +
  theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_line(color = "black", linewidth = 0.4),
        legend.position = "right") +
    annotate("text", x = 0.03 * xmax, y = max(hets_sex$het_rate) * 0.95, label = "B.")




# define boundaries
hets_ratio <- hets_sex %>%
  group_by(CHROM, win_start, win_mid) %>%
  summarise(
    n_sites = mean(n_sites, na.rm = TRUE),
    mean_female = mean(het_rate[sex == "Female"], na.rm = TRUE),
    mean_male   = mean(het_rate[sex == "Male"],   na.rm = TRUE),
    ratio = mean_female / mean_male,
    # Perform t-test: compare het_rate between females and males in this window
    pval = tryCatch(
      t.test(het_rate[sex == "Female"], het_rate[sex == "Male"])$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    padj = p.adjust(pval, method = "BH")  # FDR correction
  ) %>%
  select(CHROM, win_start, win_mid, n_sites, ratio, pval, padj)


# Combine and derive a unified significance flag
combined_ratio <- hets_ratio %>%
  mutate(
    sig_flag = dplyr::case_when(
      !is.na(padj) ~ padj < 0.001,
      TRUE         ~ FALSE
    ),
    sig_label = ifelse(sig_flag, "Significant", "Not significant")
  )

combined_ratio <- combined_ratio %>%
  arrange(win_mid) %>%
  mutate(group_id = cumsum(c(TRUE, diff(as.integer(sig_flag)) != 0)))




# 1) Ensure ordered, compute run ids for contiguous sig_flag blocks
cr <- combined_ratio %>%
  arrange(win_mid) %>%
  mutate(
    run_id = cumsum(sig_flag != dplyr::lag(sig_flag, default = first(sig_flag))) + 1
  )

# 2) Infer window size (from mid-start); used to estimate block end from starts
win_size <- 2 * stats::median(cr$win_mid - cr$win_start, na.rm = TRUE)

# 3) Summarize blocks
blocks <- cr %>%
  group_by(run_id, sig_flag) %>%
  summarise(
    block_start = min(win_start, na.rm = TRUE),
    block_end   = max(win_start, na.rm = TRUE) + win_size,
    mid         = (min(win_mid, na.rm = TRUE) + max(win_mid, na.rm = TRUE)) / 2,
    n_windows   = dplyr::n(),
    .groups = "drop"
  ) %>%
  arrange(block_start) %>%
mutate(
  block_idx = block_start / 1e6
)



# Boundaries between consecutive blocks
boundaries <- blocks %>% slice(-1L) %>% transmute(xintercept = block_start)


# Split data: multi-point runs (lines) vs singletons (segments)
line_data <- cr %>% group_by(run_id) %>% filter(dplyr::n() >= 2) %>% ungroup()

singleton_segments <- blocks %>%
  filter(n_windows == 1) %>%
  left_join(
    cr %>% group_by(run_id) %>% summarise(y = first(ratio), .groups = "drop"),
    by = "run_id"
  )

# y position for block labels
y_max   <- max(cr$ratio, na.rm = TRUE)
y_label <- y_max * 1.03

# 0) If you don't already have win_end, compute it once
win_size <- 2 * stats::median(combined_ratio$win_mid - combined_ratio$win_start, na.rm = TRUE)
cr <- combined_ratio %>%
  arrange(win_start) %>%
  mutate(win_end = win_start + win_size)

# (blocks, boundaries, y_label built as before)
p3 <- ggplot() +
  # draw each window as a horizontal segment spanning its full width
  geom_segment(
    data = cr,
    aes(x = win_start, xend = win_end, y = ratio, yend = ratio, color = sig_flag),
    linewidth = 0.9, alpha = 0.85, lineend = "butt"
  ) +
  # block boundaries
  geom_vline(
    data = boundaries,
    aes(xintercept = xintercept),
    color = "black", linewidth = 0.75, lty="11"
  ) +
  geom_vline(
    aes(xintercept = 1100000),
    color = "black", linewidth = 0.75, lty="11"
  ) +
  # block labels (use your preferred label column: block_idx or block_label)
  geom_text(
    data = blocks,
    aes(x = block_start, y = y_label, label = block_idx),  # or block_label
    inherit.aes = FALSE, color = "grey40",
    size = 2, fontface = "bold",
    angle = 90, hjust = 1, vjust = -0.2
  ) +
  labs( x = NULL,
    y = "H F:M",
    color = "Significance"
  ) +
  scale_color_manual(values = c(`FALSE` = "grey70", `TRUE` = "black"),
                     labels = c(`FALSE` = "Not significant", `TRUE` = "Significant")) +
  scale_x_continuous(limits = c(0, xmax), expand = c(0, 0)) +
  expand_limits(y = y_label * 1.02) +
  theme_bw(base_size = 15) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_line(color = "black", linewidth = 0.4),
    panel.grid.minor = element_blank(),
  ) +
  annotate("text", x = 0.03 * xmax, y = max(cr$ratio) * 0.95, label = "C.")





# collect snp density data

merged_path_snp <- "/xdisk/mcnew/dannyjackson/sulidae/analyses/neosex_inference/snpdensity_chr29/merged.snpden.tsv"
sex_path_snp    <- "/xdisk/mcnew/dannyjackson/sulidae/referencelists/sexID.csv"

# ---- read data ----
wide_snp <- read_tsv(merged_path_snp, show_col_types = FALSE)    # columns: CHROM, BIN_START, BFBO501, BFBO502, ...
sex_snp <- read_csv(
  "/xdisk/mcnew/dannyjackson/sulidae/referencelists/sexID.csv",
  col_names = c("individual", "sex"),
  show_col_types = FALSE
)

# standardize sex table column names if needed
# (tries to find a column that matches the individual IDs, then rename to 'individual')
if (!"individual" %in% names(sex_snp)) {
  cand <- intersect(names(sex_snp), c("individual","ID","id","sample","Sample","SampleID","sample_id","sampleID"))
  if (length(cand) == 0) stop("Could not find an ID column in sex table; expected something like 'individual' or 'SampleID'.")
  sex_snp <- sex_snp %>% rename(individual = !!sym(cand[1]))
}
if (!"sex" %in% names(sex_snp)) {
  cand <- intersect(names(sex_snp), c("sex","Sex","SEX"))
  if (length(cand) == 0) stop("Could not find a 'sex' column in sex table.")
  sex_snp <- sex_snp %>% rename(sex = !!sym(cand[1]))
}

# Ensure ID formatting matches column names in merged file (e.g., drop any suffixes in sex file if present)
# Example: convert any ".idxstats" -> "" in sex IDs, trim spaces
sex_snp <- sex_snp %>%
  mutate(
    individual = str_replace(individual, "\\.idxstats$", ""),
    individual = str_trim(individual)
  )

# ---- pivot to long format ----
long_snp <- wide_snp %>%
  pivot_longer(
    cols = -c(CHROM, BIN_START),
    names_to = "individual",
    values_to = "SNP_COUNT"
  ) %>%
  # if your merged used "NA" strings or blanks, coerce to numeric
  mutate(SNP_COUNT = suppressWarnings(as.numeric(SNP_COUNT)))

# ---- join sex info ----
dat_snp <- long_snp %>%
  left_join(sex_snp, by = "individual")

# Optional: set Unknown if missing
dat_snp <- dat_snp %>%
  mutate(sex = ifelse(is.na(sex) | sex == "", "Unknown", sex))

# ---- plotting ----
p4 <- ggplot(dat_snp, aes(x = BIN_START, y = SNP_COUNT,
                     group = individual, color = sex)) +
  geom_line(alpha = 0.7, linewidth = 0.4) +
  scale_color_manual(values = c(Female = "#E26D5A", Male = "#4F7CAC", Unknown = "grey50")) +
  labs( x = NULL, y = "SNPs per 100kb", color = "Sex") +
    scale_x_continuous(limits = c(0, xmax), expand = c(0, 0)) +
  theme_bw(base_size = 15) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_line(color = "black", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    legend.position = "null"
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 1.5, alpha = 1))) +
  annotate("text", x = 0.03 * xmax, y = max(dat_snp$SNP_COUNT) * 0.95, label = "A.")




#######################
# Genes of interest 
#######################


# --- Read your inputs (as you showed) ---
keratin.PC <- read_tsv("../genes_chr29/keratingenes.PC.tsv", col_names = FALSE)
colnames(keratin.PC) <- c("CHROM","gene_start","gene_end","ID","description")

isops.PC <- read_tsv("../genes_chr29/isoprenoidgenes.PC.tsv", col_names = FALSE)
colnames(isops.PC) <- c("CHROM","gene_start","gene_end","ID","description")

# --- Add source, compute midpoint, and build a combined df ---
keratin2.PC <- keratin.PC %>%
  mutate(source = "Keratin",
         mid = (gene_start + gene_end)/2)

isops2.PC <- isops.PC %>%
  mutate(source = "Isoprenoid",
         mid = (gene_start + gene_end)/2)

genesOI.PC <- bind_rows(keratin2.PC, isops2.PC) %>%
  # Create a unique row label per *source* + description so rows don't collide
  mutate(desc_row = paste(source, description, sep = ": "))

# Optional: if you're only plotting a specific chromosome/range, filter here
# genesOI <- genesOI %>% filter(CHROM == "CM062595.1", mid >= 0, mid <= xmax)

# --- Order rows: group by source, then by median midpoint (nice tidy order) ---
order_levels <- genesOI.PC %>%
  group_by(desc_row, source) %>%
  summarise(med_mid = median(mid), .groups = "drop") %>%
  arrange(source, med_mid) %>%
  pull(desc_row)

genesOI.PC <- genesOI.PC %>%
  mutate(desc_row = factor(desc_row, levels = order_levels))

# combine p5 and p2 

# Gene track PC

genesPC <- genesPC %>%
  arrange(start) %>%
  mutate(
    y = case_when(
      strand == "-" ~ 0.2,
      strand == "+" ~ 0.3,
      TRUE ~ NA_real_     # optional: handles missing strand values
    ))


genes_p5.PC <- bind_rows(keratin2.PC, isops2.PC) %>%
  mutate(y = if_else(source == "Keratin", 0.1, 0.1))


p2.PC <- ggplot() +
  # --- Gene rectangles (from p2) ---
  geom_rect(
    data = genesPC,  fill = "grey70",
    aes(xmin = start, xmax = end, ymin = y - 0.05, ymax = y + 0.05),
    color = "black", alpha = 1, linewidth = 0.1
  ) +

  # --- Points (from p5) ---
  geom_point(
    data = genes_p5.PC,
    aes(x = mid, y = y, color = source),
    size = 2, alpha = 0.9
  ) +
  scale_color_manual(
    values = c(Keratin = "#E377B2", Isoprenoid = "#41AB5D"),
    name = "Gene set"
  ) +

  # --- Shared axes and limits ---
  scale_x_continuous(limits = c(0, xmax), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.05, 0.4),
    breaks = c(0.2, 0.33),
    labels = c("-", "+")
  ) +

  # --- Labels and theme ---
  labs(
    y = NULL,
    x = NULL
  ) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    axis.text.x = element_blank()
  ) +
  annotate("text", x = 0.03 * xmax, y = max(genesPC$y) * 0.95, label = "D.")


# Gene Track MB

# --- Read your inputs (as you showed) ---
keratin.MB <- read_tsv("../genes_chr29/keratingenes.MB.tsv", col_names = FALSE)
colnames(keratin.MB) <- c("CHROM","gene_start","gene_end","ID","description")

isops.MB <- read_tsv("../genes_chr29/isoprenoidgenes.MB.tsv", col_names = FALSE)
colnames(isops.MB) <- c("CHROM","gene_start","gene_end","ID","description")

# --- Add source, compute midpoint, and build a combined df ---
keratin2.MB <- keratin.MB %>%
  mutate(source = "Keratin",
         mid = (gene_start + gene_end)/2)

isops2.MB <- isops.MB %>%
  mutate(source = "Isoprenoid",
         mid = (gene_start + gene_end)/2)

genesOI.MB <- bind_rows(keratin2.MB, isops2.MB) %>%
  # Create a unique row label per *source* + description so rows don't collide
  mutate(desc_row = paste(source, description, sep = ": "))

# Optional: if you're only plotting a specific chromosome/range, filter here
# genesOI <- genesOI %>% filter(CHROM == "CM062595.1", mid >= 0, mid <= xmax)

# --- Order rows: group by source, then by median midpoint (nice tidy order) ---
order_levels <- genesOI.MB %>%
  group_by(desc_row, source) %>%
  summarise(med_mid = median(mid), .groups = "drop") %>%
  arrange(source, med_mid) %>%
  pull(desc_row)

genesOI.MB <- genesOI.MB %>%
  mutate(desc_row = factor(desc_row, levels = order_levels))

# combine p5 and p2 

# Gene track PC

genesMB <- genesMB %>%
  arrange(start) %>%
  mutate(
    y = case_when(
      strand == "-" ~ 0.2,
      strand == "+" ~ 0.3,
      TRUE ~ NA_real_     # optional: handles missing strand values
    ))



genes_p5.MB <- bind_rows(keratin2.MB, isops2.MB) %>%
  mutate(y = if_else(source == "Keratin", 0.1, 0.1))

p2.MB <- ggplot() +
  # --- Gene rectangles (from p2) ---
  geom_rect(
    data = genesMB,  fill = "grey70",
    aes(xmin = start, xmax = end, ymin = y - 0.05, ymax = y + 0.05),
    color = "black", alpha = 1, linewidth = 0.1
  ) +

  # --- Points (from p5) ---
  geom_point(
    data = genes_p5.MB,
    aes(x = mid, y = y, color = source),
    size = 2, alpha = 0.9
  ) +
  scale_color_manual(
    values = c(Keratin = "#E377B2", Isoprenoid = "#41AB5D"),
    name = "Gene set"
  ) +

  # --- Shared axes and limits ---
  scale_x_continuous(limits = c(0, xmax), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.05, 0.4),
    breaks = c(0.2, 0.33),
    labels = c("-", "+")
  ) +

  # --- Labels and theme ---
  labs(
    y = "Strand",
    x = NULL
  ) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    axis.text.x = element_blank()
  ) +
  annotate("text", x = 0.03 * xmax, y = max(genesMB$y) * 0.95, label = "E.")

# RE Track

p2.RE <- ggplot(genesRE, aes(x = (start + end) / 2)) +
  stat_density(
    aes(y = 0.25, fill = after_stat(density)),
    geom = "tile",
    position = "identity",
    height = 0.1,
    n = 512,
    adjust = 1
  ) +
  scale_fill_gradient(low = "white", high = "black", name = "Gene density") +
  scale_x_continuous(limits = c(0, xmax), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.15, 0.35)) +
  labs(y = NULL, x = "Position (bp)") +
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    axis.text.y = element_blank()
  ) +
  annotate("text", x = 0.03 * xmax, y = 0.25, label = "F.", colour = "white")


library(grid)  # for unit()

# 1) Hide legends in each panel
plist <- list(p4, p1, p3, p2.PC, p2.MB, p2.RE)
plist_noleg <- lapply(
  plist,
  \(p) p + theme(legend.position = "none",
                 plot.margin = margin(5.5, 5.5, 2, 5.5))
)



# 2) Extract the legends you actually want to show
#    (example: one “species” legend from p1, one “Gene set” legend from p5)
leg_p1 <- get_legend(
  p1 + theme(
    legend.position = "bottom",
    legend.direction = "vertical",
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10),
    legend.box.margin = margin(t = 4, b = 4)
  )
)

leg_p2.PC <- get_legend(
  p2.PC + theme(
    legend.position = "bottom",
    legend.direction = "vertical",
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10),
    legend.box.margin = margin(t = 4, b = 4)
  )
)


leg_p2.MB <- get_legend(
  p2.MB + theme(
    legend.position = "bottom",
    legend.direction = "vertical",
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10),
    legend.box.margin = margin(t = 4, b = 4)
  )
)


leg_p2.RE <- get_legend(
  p2.RE + theme(
    legend.position = "bottom",
    legend.direction = "vertical",
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10),
    legend.box.margin = margin(t = 4, b = 4)
  )
)

leg_p3 <- get_legend(
  p3 + theme(
    legend.position = "bottom",
    legend.direction = "vertical",
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10),
    legend.box.margin = margin(t = 4, b = 4)
  )
)

leg_p4 <- get_legend(
  p4 + theme(
    legend.position = "bottom",
    legend.direction = "vertical",
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10),
    legend.box.margin = margin(t = 4, b = 4)
  )
)

# 3) Stack plots (no legends)
stacked <- plot_grid(
  plotlist = plist_noleg,
  ncol = 1, align = "v",
  rel_heights = c(0.75, 0.75, 0.75, 0.3, 0.3, 0.45)
)

# 4) Put the legends in a single row (or stack them if you prefer)
legend_row <- plot_grid(leg_p4, leg_p1, leg_p3, leg_p2.PC, leg_p2.MB, leg_p2.RE, ncol = 5)

# 5) Final combine
combined <- plot_grid(stacked, legend_row, ncol = 1, rel_heights = c(1, 0.12))

ggsave("hets_with_genes.p1234.pdf", combined, width = 12, height = 6, dpi = 300)

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(utils)
})

args  <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript pca_admix_plot.R <GROUP> <BASE_DIR>\n", call. = FALSE)
}
GROUP   <- args[1]
BASEDIR <- args[2]

# Paths
grp_dir   <- file.path(BASEDIR, GROUP)
cov_file  <- file.path(grp_dir, sprintf("genolike_pruned_%s.cov", GROUP))
# Labels live in a referencelists dir parallel to BASEDIR (per your example)
lab_file  <- "/xdisk/mcnew/dannyjackson/sulidae/referencelists/sample_species.%s.txt"
lab_file  <- sprintf(lab_file, GROUP)

if (!file.exists(cov_file)) stop("Missing covariance file: ", cov_file)
if (!file.exists(lab_file)) stop("Missing labels file: ", lab_file)

message("Reading: ", cov_file)
C <- as.matrix(read.table(cov_file, check.names = FALSE))

message("Reading: ", lab_file)
labs <- read.table(lab_file, header = TRUE, stringsAsFactors = FALSE)
# Expect columns named 'sample' and 'species'
if (!all(c("sample","species") %in% tolower(names(labs)))) {
  stop("Labels file must have columns named 'sample' and 'species' (case-insensitive).")
}
names(labs) <- tolower(names(labs))
labs$Sample  <- factor(labs$sample)
labs$Species <- factor(labs$species)

# Basic sanity check
if (nrow(labs) != nrow(C)) {
  stop(sprintf("Mismatch: %d labels vs %d rows in covariance matrix.", nrow(labs), nrow(C)))
}

# Eigen
e <- eigen(C)
PC1.PV <- round(100 * e$values[1] / sum(e$values), 2)
PC2.PV <- round(100 * e$values[2] / sum(e$values), 2)

scores <- data.frame(
  PC1 = e$vectors[, 1],
  PC2 = e$vectors[, 2],
  Sample  = labs$Sample,
  Species = labs$Species
)

# --- Make numeric labels like BFBO501 -> 1, PEBO601 -> 1, RFBO203 -> 3 ---
num_from_id <- function(x) {
  s <- sub(".*?(\\d+)$", "\\1", as.character(x))   # grab trailing digits
  n <- suppressWarnings(as.integer(s))
  lab <- n %% 100                                  # 501->1, 601->1, 203->3
  lab[lab == 0 & !is.na(lab)] <- 100               # optional: 600->100 instead of 0
  lab
}

scores$LabelNum <- num_from_id(scores$Sample)


# Shared helpers
# Make panel truly square and visually balanced
square_panel <- function(p, df = scores) {
  rng <- range(c(df$PC1, df$PC2), finite = TRUE)
  p + coord_fixed(ratio = 1, xlim = rng, ylim = rng)
}

# Species color map (extendable; unknowns -> grey70)
species_cols <- c(
  BFBO = "#4FE98C",
  PEBO = "#2894FF",
  RFBO = "#F3447D",
  BRBO = "#000000",
  MABO = "#E8CA3D",
  NABO = "#B32904"
)
# Add any species present but not in map as grey
missing_sp <- setdiff(levels(scores$Species), names(species_cols))
if (length(missing_sp)) {
  add <- stats::setNames(rep("grey70", length(missing_sp)), missing_sp)
  species_cols <- c(species_cols, add)
}

# ===== PCA: species-colored =====
p_species <- ggplot(scores, aes(x = PC1, y = PC2, fill = Species)) +
  geom_point(shape = 21, size = 10, alpha = 0.9, color = "white", stroke = 0.3) +
  scale_fill_manual(values = species_cols, drop = FALSE) +
  labs(
    x = paste0("PC1 (", PC1.PV, "%)"),
    y = paste0("PC2 (", PC2.PV, "%)")
  ) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  )

p_species <- square_panel(p_species, scores)

out_species <- file.path(grp_dir, sprintf("pca_pruned_%s.species.pdf", GROUP))
ggsave(out_species, p_species, width = 5.5, height = 5.5, device = cairo_pdf)
message("Wrote: ", out_species)

# ===== PCA: sample-colored (big legends can be unwieldy; keep it) =====
p_sample <- ggplot(scores, aes(x = PC1, y = PC2, fill = Sample)) +
  geom_point(shape = 21, size = 15, alpha = 0.9, color = "white", stroke = 0.3) +
  labs(
    x = paste0("PC1 (", PC1.PV, "%)"),
    y = paste0("PC2 (", PC2.PV, "%)"),
    fill = "Sample"
  ) +
  theme_bw(base_size = 21) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    legend.key.height = grid::unit(0.35, "lines"),
    legend.text = element_text(size = 7)
  )

p_sample <- square_panel(p_sample, scores)

out_sample <- file.path(grp_dir, sprintf("pca_pruned_%s.sample.pdf", GROUP))
ggsave(out_sample, p_sample, width = 5.5, height = 5.5, device = cairo_pdf)
message("Wrote: ", out_sample)

# ===== ADMIXTURE barplot =====
# find a *.Q file for this group
q_candidates <- Sys.glob(file.path(grp_dir, sprintf("genolike_pruned_%s.admix*.Q", GROUP)))
if (length(q_candidates) == 0L) {
  message("No admixture Q file found for group ", GROUP, "; skipping admixture plot.")
  quit(status = 0)
}
q_file <- q_candidates[1]
message("Reading admixture: ", q_file)

Q <- read.table(q_file, header = FALSE, stringsAsFactors = FALSE)
if (nrow(Q) != nrow(labs)) {
  warning(sprintf("Admixture individuals (%d) != labels (%d). Proceeding, but order may differ.",
                  nrow(Q), nrow(labs)))
}
# If you want bars in the same order as labels file:
ord <- seq_len(nrow(Q))
# (If your Q is in label order, this is fine. If not, match here.)

out_admix <- file.path(grp_dir, sprintf("admix_pruned_%s.pdf", GROUP))
pdf(file = out_admix, width = 8, height = 4, family = "Helvetica")
par(mar = c(4, 4, 1, 1) + 0.1)
barplot(t(as.matrix(Q[ord, ])),
        col = grDevices::rainbow(ncol(Q)),
        xlab = "Individual #", ylab = "Ancestry",
        border = NA, axes = TRUE)
dev.off()
message("Wrote: ", out_admix)

message("Done: ", GROUP)

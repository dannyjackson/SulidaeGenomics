#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(stringr)
  library(rrvgo)
  library(org.Gg.eg.db)
})

## --------------------------------------------------------------------
## 1. Load and harmonize data
## --------------------------------------------------------------------

# NOTE: switch to read_tsv() if your files are tab-delimited
brbo <- readr::read_csv("BRBO.fixedsites_in_genes.GO.Gallus.newhead.fdr.2.txt")
mabo <- readr::read_csv("MABO.fixedsites_in_genes.GO.Gallus.newhead.fdr.2.txt")

# Extract GO IDs from GO_Term strings
extract_id <- function(x) stringr::str_extract(x, "GO:\\d+")
brbo$GO_ID <- extract_id(brbo$GO_Term)
mabo$GO_ID <- extract_id(mabo$GO_Term)

# Keep GO terms with at least one gene
brbo <- brbo %>% dplyr::filter(Empirical > 0)
mabo <- mabo %>% dplyr::filter(Empirical > 0)

# Shared GO IDs
shared_ids <- intersect(brbo$GO_ID, mabo$GO_ID)

brbo_filt <- brbo %>%
  dplyr::filter(GO_ID %in% shared_ids) %>%
  dplyr::mutate(Species = "BRBO")

mabo_filt <- mabo %>%
  dplyr::filter(GO_ID %in% shared_ids) %>%
  dplyr::mutate(Species = "MABO")

# Make columns numeric and consistent
brbo_filt <- brbo_filt %>%
  dplyr::mutate(
    P.value = as.numeric(P.value),
    # handle things like "> 100" or spaces
    foldEnrichment = stringr::str_replace_all(foldEnrichment, ">", ""),
    foldEnrichment = stringr::str_replace_all(foldEnrichment, "\\s+", ""),
    foldEnrichment = as.numeric(foldEnrichment),
    qvals = as.numeric(qvals)
  )

mabo_filt <- mabo_filt %>%
  dplyr::mutate(
    P.value        = as.numeric(P.value),
    foldEnrichment = as.numeric(foldEnrichment),
    qvals          = as.numeric(qvals)
  )

# Combine into one data frame
df <- dplyr::bind_rows(brbo_filt, mabo_filt) %>%
  dplyr::mutate(
    FDR = as.numeric(qvals),
    # avoid FDR == 0 causing -Inf
    FDR = dplyr::if_else(FDR == 0, 1e-320, FDR)
  )

## --------------------------------------------------------------------
## 2. Shared GO-term dotplot (no clustering)
## --------------------------------------------------------------------

p_simple <- ggplot(df, aes(x = Species,
                           y = GO_Term,
                           size = foldEnrichment,
                           color = -log10(FDR))) +
  geom_point(alpha = 0.85) +
  scale_color_gradient(low = "yellow", high = "red") +
  scale_size(range = c(1, 4)) +
  labs(x = "", y = "",
       color = "-log10(FDR)",
       size  = "Fold enrichment") +
  theme_bw(base_size = 8) +
  theme(
    axis.text.y        = element_text(size = 8),
    axis.text.x        = element_text(size = 10),
    panel.grid.major.y = element_blank()
  )

ggsave("GO_shared_dotplot.pdf", p_simple, width = 8, height = 12)

## --------------------------------------------------------------------
## 3. Cluster GO terms by function using rrvgo
## --------------------------------------------------------------------

# Build a score per GO_ID (max -log10(FDR) across species)
scores_df <- df %>%
  dplyr::group_by(GO_ID) %>%
  dplyr::summarise(score = max(-log10(FDR), na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(score = dplyr::if_else(is.infinite(score), NA_real_, score)) %>%
  dplyr::filter(!is.na(score))

scores <- scores_df$score
names(scores) <- scores_df$GO_ID

# Only run rrvgo if we have at least 2 GO terms with scores

simMatrix <- rrvgo::calculateSimMatrix(
  names(scores),
  orgdb  = "org.Gg.eg.db",
  ont    = "BP",
  method = "Rel"
)

reducedTerms <- rrvgo::reduceSimMatrix(
  simMatrix,
  scores,
  orgdb     = "org.Gg.eg.db",
  threshold = 0.4
)

# Map each original GO ID (term) to its parent cluster name
cluster_map <- reducedTerms %>%
  dplyr::transmute(
    GO_ID       = go,        # original GO ID
    ClusterName = parentTerm   # functional cluster label
  )

# Join cluster info back to df
df_plot <- df %>%
  dplyr::left_join(cluster_map, by = "GO_ID")

## ------------------------------------------------------------------
## 4. Order GO terms by cluster + significance
## ------------------------------------------------------------------

order_df <- df_plot %>%
  dplyr::select(GO_ID, GO_Term, ClusterName) %>%
  dplyr::distinct() %>%
  dplyr::left_join(scores_df, by = "GO_ID") %>%
  dplyr::arrange(ClusterName, dplyr::desc(score))

term_levels <- order_df$GO_Term

df_plot <- df_plot %>%
  dplyr::mutate(
    GO_Term = factor(GO_Term, levels = term_levels)
  )

## ------------------------------------------------------------------
## 5. Cluster membership bar plot (all GO terms on x-axis)
## ------------------------------------------------------------------

# One row per GO term: which cluster is it in?
bar_df <- df_plot %>%
  dplyr::select(GO_Term, ClusterName) %>%
  dplyr::distinct() %>%
  dplyr::mutate(GO_Term = factor(GO_Term, levels = term_levels),
                value   = 1)

# This gives a single bar per GO term, colored by cluster
p_bar <- ggplot(bar_df, aes(x = GO_Term, y = value, fill = ClusterName)) +
  geom_col() +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "GO term", y = "Cluster membership") +
  theme_bw(base_size = 8) +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid   = element_blank()
  )

ggsave("GO_cluster_barplot.pdf", p_bar, width = 10, height = 4)


p_simple <- ggplot(df, aes(x = Species,
                           y = GO_Term,
                           size = foldEnrichment,
                           color = -log10(FDR))) +
  geom_point(alpha = 0.85) +
  scale_color_gradient(low = "yellow", high = "red") +
  scale_size(range = c(1, 4)) +
  labs(x = "", y = "",
       color = "-log10(FDR)",
       size  = "Fold enrichment") +
  theme_bw(base_size = 8) +
  theme(
    axis.text.y        = element_text(size = 8),
    axis.text.x        = element_text(size = 10),
    panel.grid.major.y = element_blank()
  )

ggsave("GO_shared_dotplot.pdf", p_simple, width = 8, height = 12)




# One row per GO term: which cluster is it in?
df_cluster <- df_plot %>%
  dplyr::select(GO_Term, ClusterName) %>%
  dplyr::distinct() %>%
  dplyr::mutate(SpeciesPanel = "Cluster")

# Species-level data (BRBO / MABO)
df_species <- df_plot %>%
  dplyr::mutate(SpeciesPanel = Species)
  
## 5. Combined plot: cluster point + dotplot ----
p_combined <- ggplot() +
  # Cluster-colored point (one per GO term)
  geom_point(
    data = df_cluster,
    aes(x = SpeciesPanel, y = GO_Term, color = ClusterName),
    size = 2
  ) +
  # Species × GO dotplot, colored by FDR, sized by foldEnrichment
  geom_point(
    data = df_species,
    aes(x = SpeciesPanel, y = GO_Term,
        size = foldEnrichment,
        fill = -log10(FDR)),
    shape = 21,
    color = "black",
    alpha = 0.85
  ) +
  scale_x_discrete(limits = c("Cluster", "BRBO", "MABO"), name = "") +
  scale_fill_gradient(low = "yellow", high = "red", name = "-log10(FDR)") +
  scale_color_discrete(name = "Cluster") +
  scale_size(range = c(1, 4), name = "Fold enrichment") +
  labs(y = "") +
  theme_bw(base_size = 8) +
  theme(
    axis.text.y        = element_text(size = 7),
    axis.text.x        = element_text(size = 9),
    panel.grid.major.y = element_blank()
  )

ggsave("GO_shared_dotplot_with_clusters.pdf", p_combined, width = 8, height = 12)


# df_cluster has one row per GO term:
# GO_Term | ClusterName | SpeciesPanel == "Cluster"
cluster_blocks <- df_cluster %>%
  dplyr::mutate(y_numeric = as.numeric(GO_Term)) %>%
  dplyr::group_by(ClusterName) %>%
  dplyr::summarise(
    y_min = min(y_numeric) - 0.45,
    y_max = max(y_numeric) + 0.45,
    .groups = "drop"
  )

p_combined <- ggplot() +
  geom_rect(
    data = cluster_blocks,
    aes(xmin = 0, xmax = 0.3, ymin = y_min, ymax = y_max, fill = ClusterName),
    color = NA,
    alpha = 0.9
  ) +
  geom_point(
    data = df_species,
    aes(
      x = factor(SpeciesPanel, levels = c("Cluster", "BRBO", "MABO")),
      y = GO_Term,
      size = foldEnrichment,
      color = -log10(FDR)
    ),
    alpha = 0.85
  ) +
  scale_x_discrete(limits = c("Cluster", "BRBO", "MABO"), name = "") +
  scale_fill_discrete(name = "Cluster") +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(FDR)") +
  scale_size(range = c(1, 4), name = "Fold enrichment") +
  labs(y = "") +
  theme_bw(base_size = 8) +
  theme(
    axis.text.y        = element_text(size = 7),
    axis.text.x        = element_text(size = 9),
    panel.grid.major.y = element_blank(),

    # 🔥 Legend shrink controls:
    legend.key.size    = unit(0.4, "lines"),
    legend.text        = element_text(size = 6),
    legend.title       = element_text(size = 6),
    legend.spacing.y   = unit(0.1, "lines"),
    legend.spacing.x   = unit(0.2, "lines"),

    # optional: move to bottom to avoid shrinking plotting space
    legend.position    = "bottom"
  )

ggsave("GO_shared_dotplot_with_merged_cluster_blocks.pdf", p_combined, width = 6, height = 12)

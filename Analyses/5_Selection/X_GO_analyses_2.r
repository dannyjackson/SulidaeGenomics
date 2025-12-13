#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(stringr)
  library(rrvgo)
  library(GOSemSim)
  library(org.Gg.eg.db)
})

## --------------------------------------------------------------------
## 1. Load and harmonize data
## --------------------------------------------------------------------

# Change read_csv() to read_tsv() if your files are tab-delimited
brbo <- readr::read_csv("BRBO.fixedsites_in_genes.GO.Gallus.newhead.fdr.2.txt")
mabo <- readr::read_csv("MABO.fixedsites_in_genes.GO.Gallus.newhead.fdr.2.txt")

# Extract GO IDs from GO_Term strings
extract_id <- function(x) stringr::str_extract(x, "GO:\\d+")
brbo$GO_ID <- extract_id(brbo$GO_Term)
mabo$GO_ID <- extract_id(mabo$GO_Term)

# Keep GO terms with at least one gene
brbo <- brbo %>% filter(Empirical > 0)
mabo <- mabo %>% filter(Empirical > 0)

# Shared GO IDs
shared_ids <- intersect(brbo$GO_ID, mabo$GO_ID)

brbo_filt <- brbo %>%
  filter(GO_ID %in% shared_ids) %>%
  mutate(Species = "BRBO")

mabo_filt <- mabo %>%
  filter(GO_ID %in% shared_ids) %>%
  mutate(Species = "MABO")

# Make columns numeric and consistent
brbo_filt <- brbo_filt %>%
  mutate(
    P.value = as.numeric(P.value),
    # handle things like "> 100" or spaces
    foldEnrichment = stringr::str_replace_all(foldEnrichment, ">", ""),
    foldEnrichment = stringr::str_replace_all(foldEnrichment, "\\s+", ""),
    foldEnrichment = as.numeric(foldEnrichment),
    qvals = as.numeric(qvals)
  )

mabo_filt <- mabo_filt %>%
  mutate(
    P.value = as.numeric(P.value),
    foldEnrichment = as.numeric(foldEnrichment),
    qvals = as.numeric(qvals)
  )

# Combine into one data frame
df <- dplyr::bind_rows(brbo_filt, mabo_filt) %>%
  # use qvals as FDR
  mutate(FDR = as.numeric(qvals)) %>%
  # avoid FDR == 0 causing -Inf
  mutate(FDR = ifelse(FDR == 0, 1e-320, FDR))

## --------------------------------------------------------------------
## 2. Simple shared GO-term dotplot (no clustering)
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
  group_by(GO_ID) %>%
  summarise(score = max(-log10(FDR), na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(score))

scores <- scores_df$score
names(scores) <- scores_df$GO_ID

# Semantic similarity matrix for GO IDs
simMatrix <- rrvgo::calculateSimMatrix(
  names(scores),
  orgdb  = "org.Gg.eg.db",
  ont    = "BP",
  method = "Wang"
)

# Reduce to parent terms (functional clusters)
reducedTerms <- rrvgo::reduceSimMatrix(
  simMatrix,
  scores,
  orgdb     = "org.Gg.eg.db",
  threshold = 0.7  # lower = bigger clusters, higher = more granular
)

# Map each GO_ID to its parent cluster
cluster_map <- reducedTerms %>%
  transmute(
    GO_ID       = parent,
    ClusterID   = term,
    ClusterName = parentTerm
  )

ClusterID
# Join cluster info back to df
df_plot <- df %>%
  inner_join(cluster_map, by = "GO_ID")  # drop terms not in reducedTerms

## --------------------------------------------------------------------
## 4. Order GO terms within clusters by significance for plotting
## --------------------------------------------------------------------

# Reuse scores_df (by GO_ID) and join to cluster info
order_df <- df_plot %>%
  dplyr::select(GO_ID, GO_Term, ClusterName) %>%
  dplyr::distinct() %>%
  dplyr::left_join(scores_df, by = "GO_ID") %>%
  dplyr::arrange(ClusterName, dplyr::desc(score))

# Order of GO_Term levels
term_levels <- order_df$GO_Term

# Order of cluster facets (in same order they appear in order_df)
cluster_levels <- order_df %>%
  distinct(order_df$ClusterName) %>%
  pull(ClusterName)

df_plot <- df_plot %>%
  mutate(
    GO_Term     = factor(GO_Term, levels = term_levels),
    ClusterName = factor(ClusterName, levels = cluster_levels)
  )

## --------------------------------------------------------------------
## 5. Faceted dotplot by functional cluster
## --------------------------------------------------------------------

p_clustered <- ggplot(df_plot, aes(x = Species,
                                   y = GO_Term,
                                   size = foldEnrichment,
                                   color = -log10(FDR))) +
  geom_point(alpha = 0.85) +
  scale_color_gradient(low = "yellow", high = "red") +
  scale_size(range = c(2, 12)) +
  labs(x = "", y = "",
       color = "-log10(FDR)",
       size  = "Fold enrichment") +
  facet_grid(
    rows   = vars(ClusterName),
    scales = "free_y",
    space  = "free_y"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y        = element_text(size = 8),
    axis.text.x        = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    strip.placement    = "outside",
    strip.background   = element_rect(fill = "grey90"),
    strip.text.y.left  = element_text(angle = 0, face = "bold")
  )

ggsave("GO_shared_dotplot_rrvgo_clusters.pdf", p_clustered, width = 10, height = 12)

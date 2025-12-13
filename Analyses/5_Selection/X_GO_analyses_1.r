#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(ggplot2)
library(stringr)

# ---- Load files ----
brbo <- readr::read_csv("BRBO.fixedsites_in_genes.GO.Gallus.newhead.fdr.2.txt")
mabo <- readr::read_csv("MABO.fixedsites_in_genes.GO.Gallus.newhead.fdr.2.txt")

# ---- Extract the GO IDs from GO_Term column ----
extract_id <- function(x) str_extract(x, "GO:\\d+")
brbo$GO_ID <- extract_id(brbo$GO_Term)
mabo$GO_ID <- extract_id(mabo$GO_Term)

# ---- Keep GO terms with at least one gene (Empirical > 0) ----
brbo <- brbo %>% filter(Empirical > 0)
mabo <- mabo %>% filter(Empirical > 0)

# ---- Shared GO terms ----
shared <- intersect(brbo$GO_ID, mabo$GO_ID)

brbo_filt <- brbo %>% filter(GO_ID %in% shared) %>% mutate(Species = "BRBO")
mabo_filt <- mabo %>% filter(GO_ID %in% shared) %>% mutate(Species = "MABO")


# Make BRBO numeric where needed
brbo_filt <- brbo_filt %>%
  mutate(
    # handle things like "3.40E-03" safely
    P.value = as.numeric(P.value),

    # deal with possible "> 100" etc
    foldEnrichment = if_else(
      str_detect(foldEnrichment, ">"),
      "100",
      foldEnrichment
    ),
    foldEnrichment = as.numeric(foldEnrichment),

    # qvals currently character, make numeric
    qvals = as.numeric(qvals)
  )

# Make MABO columns the same types
mabo_filt <- mabo_filt %>%
  mutate(
    P.value        = as.numeric(P.value),        # already double, but harmless
    foldEnrichment = as.numeric(foldEnrichment),
    qvals          = as.numeric(qvals)
  )

# ---- Combine ----
df <- bind_rows(brbo_filt, mabo_filt)

# Use qvals as FDR
df <- df %>%
  mutate(FDR = qvals)


# ---- FDR correction (per species) ----
df <- df %>%
  group_by(Species) %>%
  mutate(FDR = p.adjust(P.value, method = "fdr")) %>%
  ungroup()

# ---- Convert foldEnrichment to numeric if has "> 100" ----
df$foldEnrichment <- as.numeric(str_replace(df$foldEnrichment, "> 100", "100"))

# ---- Plot ----
p <- ggplot(df, aes(x = Species,
                    y = GO_Term,
                    size = foldEnrichment,
                    color = -log10(FDR))) +
  geom_point(alpha = 0.85) +
  scale_color_gradient(low = "yellow", high = "red") +
  scale_size(range = c(1, 4)) +
  labs(x = "", y = "",
       color = "-log10(FDR)",
       size = "Fold enrichment") +
  theme_bw(base_size = 8) +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 10),
    panel.grid.major.y = element_blank()
  )

ggsave("GO_shared_dotplot.pdf", p, width = 8, height = 12)


# Use rrvgo to parent terms
library(dplyr)
library(rrvgo)

library(GOSemSim)
library(org.Gg.eg.db)
library(cluster)

# extract GO IDs
go_ids <- unique(df$GO_ID)

# build ontology object
godata_obj <- godata("org.Gg.eg.db", ont = "BP")

# similarity matrix
sim_mat <- mgeneSim(go_ids, semData = godata_obj, measure = "Wang")

# hierarchical clustering
hc <- hclust(as.dist(1 - sim_mat), method = "average")

# cluster assignment (change k if you want more/less clusters)
cluster_assign <- cutree(hc, k = 6)

cluster_df <- data.frame(GO_ID = names(cluster_assign),
                         Cluster = factor(cluster_assign))




# 1) Make a named score vector for each GO term
#    You can use -log10(FDR), foldEnrichment, or whatever you prefer.
scores_df <- df %>%
  group_by(GO_Term) %>%
  summarise(score = max(-log10(FDR), na.rm = TRUE), .groups = "drop")

scores <- scores_df$score
names(scores) <- scores_df$GO_Term

# 2) Semantic similarity matrix using rrvgo helper
#    Change orgdb/ont to your organism & ontology
simMatrix <- calculateSimMatrix(
  names(scores),
  orgdb  = "org.Hs.eg.db",
  ont    = "BP",
  method = "Wang"
)

# 3) Reduce to higher-order "parent" terms (this is what treePlot() uses)
reducedTerms <- reduceSim(
  simMatrix,
  scores,
  orgdb     = "org.Gg.eg.db",
  threshold = 0.7  # tweak: lower = bigger clusters, higher = more granular
)

head(reducedTerms)
# columns include:
#   term        (original GO ID)
#   termName    (original term name)
#   parent      (representative GO ID)
#   parentTerm  (representative term name)
#   score, uniqueness, dispensability, ...


# Step 2 – Map each GO term in your df to its parent term
cluster_map <- reducedTerms %>%
  transmute(
    GO_Term     = term,
    ClusterID   = parent,
    ClusterName = parentTerm
  )

df_plot <- df %>%
  left_join(cluster_map, by = "GO_Term")

# Step 3 – Order GO terms by cluster, then within-cluster
# order terms within clusters by decreasing score (or whatever)
term_order <- reducedTerms %>%
  arrange(ClusterName = parentTerm, desc(score)) %>%
  pull(term)

# order cluster labels by where they appear along the y-axis
cluster_order <- reducedTerms %>%
  mutate(term = factor(term, levels = term_order),
         y = as.numeric(term)) %>%
  group_by(parentTerm) %>%
  summarise(mean_y = mean(y), .groups = "drop") %>%
  arrange(desc(mean_y)) %>%
  pull(parentTerm)

df_plot <- df_plot %>%
  mutate(
    GO_Term     = factor(GO_Term, levels = term_order),
    ClusterName = factor(ClusterName, levels = cluster_order)
  )

# Step 4 – Use facets as “bars” for each parent term along the y-axis
library(ggplot2)

p <- ggplot(df_plot, aes(x = Species,
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

ggsave("GO_shared_dotplot_rrvgo_clusters.pdf", p, width = 10, height = 12)

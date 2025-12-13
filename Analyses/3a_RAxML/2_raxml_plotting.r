# Make plots
library(ape)
library(phangorn)
library(ggplot2)
library(ggtree)
library(treeio)
library(patchwork)
library(scales)
library(readr)
library(stringr)
library(phytools)


# Read all trees
trees <- read.tree("allbestTree.autosomes.nwk")  # or list.files() and lapply(read.tree, ...)

# plot distance 
# Compute pairwise Robinson-Foulds distances
rf <- dist.topo(trees, method="PH85")

# Ordinate in 2D
mds <- cmdscale(rf, k=2)
mds_df <- data.frame(mds, chromosome = c(1:31))

ggplot(mds_df, aes(X1, X2, label=chromosome)) +
  geom_point() +
  geom_text(vjust=-0.5) +
  theme_minimal() +
  labs(title="Chromosome Tree Topology Distances (RF metric)")

ggsave("autosome.bestTree.distance.png")
dev.off()


# ---- 1) Load your 31 chromosome trees ----
dir_trees <- "/xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/vcfs/trees/autosomes/"
files <- Sys.glob(file.path(dir_trees, "*bestTree"))
# If needed, still exclude any:

trees <- lapply(files, read.tree)

# ---- 2) Root on 'reference' and ladderize for consistency ----
trees_rooted <- lapply(trees, function(tr) {
  if ("reference" %in% tr$tip.label) {
    tr <- root(tr, outgroup = "reference", resolve.root = TRUE)
  } else {
    warning("A tree is missing the tip 'reference' — left unrooted.")
  }
  ladderize(tr)
})

# ---- Build tip annotation (shared across trees if same taxa set) ----
tip_labels <- trees_rooted[[1]]$tip.label


# species prefix = all chars before trailing digits (e.g., BFBO from BFBO501)
species <- sub("(.*?)(\\d+)$", "\\1", tip_labels)       # strip trailing digits
species <- sub("_.*$", "", species)                     # drop anything after underscores if present
species[tip_labels == "reference"] <- "reference"

# last digit from the sample name; blank for 'reference'
last_digit <- ifelse(tip_labels == "reference", "ref",
                     sub("^.*?(\\d)$", "\\1", tip_labels))

tip_df <- data.frame(label = tip_labels,
                     species = species,
                     digit = last_digit,
                     stringsAsFactors = FALSE)


# ---- Read fixed colors from TSV ----
cc <- read_tsv("/xdisk/mcnew/dannyjackson/sulidae/referencelists/colorcodes.tsv", col_types = "cc", progress = FALSE)
# make a named vector like c("BFBO"="#4FE98C", ...)
fixed_map <- setNames(cc$Code, toupper(cc$Species))

# ---- Build species list present in the plots ----
species_levels <- sort(unique(tip_df$species))  # tip_df$species should be like BFBO/PEBO/etc.

# ---- Start palette with NA, fill from file, then handle special/missing ----
pal <- setNames(rep(NA_character_, length(species_levels)), species_levels)

# fill known species from the file
overlap <- intersect(names(fixed_map), species_levels)
pal[overlap] <- fixed_map[overlap]

# always set reference to grey (if present)
if ("reference" %in% names(pal)) pal["reference"] <- "grey55"

# ---- derive sites per file from matching *.raxml.log ----
get_chr <- function(f) {
  b <- basename(f)
  b <- sub("\\.raxml\\.bestTree$", "", b)
  b <- sub("\\.bestTree$", "", b)
  b
}

get_sites_from_log <- function(chr, log_dir) {
  logf <- file.path(log_dir, paste0(chr, ".raxml.log"))
  if (!file.exists(logf)) return(NA_integer_)
  ln <- grep("Alignment sites[[:space:]]*/[[:space:]]*patterns", readLines(logf, warn = FALSE), value = TRUE)
  if (!length(ln)) return(NA_integer_)
  # extract the first number before the slash
  sites <- sub(".*:\\s*([0-9,]+)\\s*/.*", "\\1", ln[1])
  as.integer(gsub(",", "", sites))
}


chroms_df <- read.csv("summary.tsv", sep="\t")

chroms <- chroms_df$sites

# Compute % support from bootstrap trees, write into node.label, and collapse weak edges
add_bs_and_collapse <- function(tree_file, thr = 50) {
  tr <- read.tree(tree_file)

  # locate bootstrap tree file (adjust pattern if needed)
  bs_file <- sub("\\.raxml\\.bestTree$|\\.bestTree$", ".raxml.bootstraps", tree_file)

  if (file.exists(bs_file)) {
    bs_trees <- read.tree(bs_file)  # multiphylo list of bootstrap trees
    # support in percent for each internal node of 'tr'
    supp <- ape::prop.clades(tr, bs_trees)/2
    tr$node.label <- sprintf("%d", round(supp))  # store as character for ggtree

    # collapse branches with support < thr
    if (!is.null(tr$edge.length)) {
      low_nodes <- which(supp < thr) + Ntip(tr)  # internal node numbers
      if (length(low_nodes)) {
        idx <- which(tr$edge[,2] %in% low_nodes)
        if (length(idx)) {
          tr$edge.length[idx] <- 0
          tr <- ape::di2multi(tr, tol = 0)       # collapse to polytomies
        }
      }
    }
  } else {
    warning("No bootstrap file for: ", basename(tree_file), " — skipping BS/collapse")
  }

  tr
}

THR <- 70  # change to 70 if you prefer

# Rebuild trees with BS + collapse, then root & ladderize
trees_bs <- lapply(files, add_bs_and_collapse, thr = THR)
trees_rooted <- lapply(trees_bs, function(tr) {
  if ("reference" %in% tr$tip.label) tr <- root(tr, outgroup = "reference", resolve.root = TRUE)
  ladderize(tr)
})
#        geom_text2(aes(subset = (label != "reference"), label = digit, color = species),
                  # nudge_x = 0.002, nudge_y = 0, size = 1.5) +
# ---- rebuild plots (titles/subtitles already set earlier) ----
plots <- Map(function(tr, panel_num, subtitle_txt) {
  p <- ggtree(tr, size = 0.25) %<+% tip_df +

       geom_tiplab(                          # adds name of individual to tip of its branch
            aes(color = species, label = digit),                    
            linesize = 0.25,
            geom = "text",
            align = TRUE,
            size = 2,
            linetype = "dotted",
            show.legend = FALSE,
            offset = 0.005
            ) +
            
        geom_tippoint(aes(color = species), size = 0.5, alpha = 0.9, position = position_nudge(x = 0.0015) )


    # HOLLOW nodes for 0–70
    p <- p + geom_point2(
    aes(subset = !isTip &
                    !is.na(as.numeric(label)) &
                    as.numeric(label) < THR),
    shape = 23,            # circle w/ fill + border
    size = 0.75,
    stroke = 0.3,          # border thickness
    fill = "white",
    color = "black"
    )

    # GRADIENT-FILLED nodes for 70–100
    p <- p + geom_point2(
    aes(subset = !isTip &
                    !is.na(as.numeric(label)) &
                    as.numeric(label) >= THR,
        fill = as.numeric(label)),
    shape = 23,            # circle w/ fill + border
    size = 0.75,
    stroke = 0.3,
    color = "black"
    ) +
    scale_fill_gradient(
        name   = "Bootstrap",
        limits = c(THR, 100),
        low    = "white",
        high   = "black",
        na.value = NA
    )
  
  xr <- range(p$data$x, na.rm = TRUE)
  p <- p + coord_cartesian(xlim = c(xr[1], xr[2] * 1.08), clip = "off")

  p + labs(
      title    = if (panel_num == 34) "Chromosome Z" else paste("Chromosome", panel_num),
      subtitle = paste("Number of sites:", subtitle_txt)
    ) +
  scale_color_manual(values = pal, breaks = species_levels, name = "Species") +
  theme(
    plot.title    = element_text(size = 8, hjust = 0.5),
    plot.subtitle = element_text(size = 6, hjust = 0.5, margin = margin(b = 2)),
    plot.margin   = margin(1, 4, 4, 1)
  )

}, trees_rooted, seq_along(trees_rooted), chroms)



# ---- collect legend and save ----
# make both legends vertical (stacked)
guides_vertical <- guides(
  color = guide_legend(title = "Species", ncol = 1, byrow = TRUE),
  fill  = guide_legend(title = "Bootstrap (%)", ncol = 1, byrow = TRUE,
                       override.aes = list(shape = 21, size = 3, color = "black"))
)


ncol <- 5
nr <- ceiling(length(plots) / ncol)
nslots <- nr * ncol
nempty <- nslots - length(plots)

# pad so the legend occupies the bottom-right cell of the grid
if (nempty > 0) {
  pad <- replicate(max(nempty - 1, 0), plot_spacer(), simplify = FALSE)
  pad <- c(pad, list(guide_area()))   # legend in the very last cell
  plots_padded <- c(plots, pad)
} else {
  # if no empty cells, add a new legend-only row/col (rare for your case)
  plots_padded <- c(plots, guide_area())
  nr <- ceiling(length(plots_padded) / ncol)
}

grid <- wrap_plots(plots_padded, ncol = ncol) +
  plot_layout(guides = "collect")

# make both legends vertical + place at bottom-right of that cell

grid <- grid &
  guides(
    color = guide_legend(
      title = "Species", ncol = 1, byrow = TRUE, order = 1,
      keyheight = unit(6, "pt"), keywidth = unit(6, "pt")
    ),
    fill  = guide_legend(
      title = "Bootstrap (%)", ncol = 1, byrow = TRUE, order = 2,
      override.aes = list(shape = 21, size = 2.5, color = "black"),
      keyheight = unit(6, "pt"), keywidth = unit(6, "pt")
    )
  ) &
  theme(
    # put the legend where you want it
    legend.position   = "right",
    legend.direction  = "vertical",
    legend.box        = "vertical",

    # tighten spacing BETWEEN the two legends
    legend.spacing.y  = unit(2, "pt"),

    # tighten spacing AROUND the legend box
    legend.box.spacing = unit(0, "pt"),
    legend.margin      = margin(0, 0, 0, 0),

    # tighten spacing WITHIN each legend (between keys/rows)
    legend.key.height = unit(6, "pt"),
    legend.key.width  = unit(6, "pt"),
    legend.spacing    = unit(1, "pt"),

    # text sizes
    legend.title = element_text(size = 8),
    legend.text  = element_text(size = 6),

    # cleaner background
    legend.background = element_rect(fill = NA, colour = NA)
  )

# solve versioning issue of is.waive being depreciated in newer ggplot2

if (!exists("is.waive")) {
  is.waive <- function(x) inherits(x, "waiver")
}


pdf("chromosome_trees_grid.rooted.colored_digits.legend_sites.pdf",
    width = 7, height = 20, useDingbats = FALSE)
print(grid)
dev.off()

ggsave("chromosome_trees_grid.rooted.colored_digits.legend_sites.png",
       grid, width = 7, height = 10, dpi = 300)

    




# generate consensus autosomal tree

# First, evaluate distance with and without 33

trees <- read.tree("allbestTree.autosomes.nwk")  # or list.files() and lapply(read.tree, ...)


# plot distance 
# Compute pairwise Robinson-Foulds distances
rf <- dist.topo(trees, method="PH85")

# Ordinate in 2D
mds <- cmdscale(rf, k=2)
mds_df <- data.frame(mds, chromosome = c(1:31))

ggplot(mds_df, aes(X1, X2, label=chromosome)) +
  geom_point() +
  geom_text(vjust=-0.5) +
  theme_minimal() +
  labs(title="Chromosome Tree Topology Distances (RF metric)")

ggsave("autosomes.bestTree.distance.pdf")
dev.off()




# Compute consensus network
# Make plots

net <- consensusNet(trees, prob = 0.1)

pdf(file = "autosomes.bestTree.pdf", width = 10, height = 10, useDingbats=FALSE)
plot(net, main="Consensus Network Across Chromosomes")
dev.off()

# make consensus tree

trees <- read.tree("allbestTree.autosomes.nwk")  # or list.files() and lapply(read.tree, ...)

cons <- consensus(trees, p=0.5)  # majority-rule consensus
cons_rooted <- root(cons, outgroup="reference", resolve.root=TRUE)
labs <- round(100 * prop.clades(cons_rooted, trees), 1)


# Build species color code matrix

tip_labels <- cons_rooted$tip.label


# species prefix = all chars before trailing digits (e.g., BFBO from BFBO501)
species <- sub("(.*?)(\\d+)$", "\\1", tip_labels)       # strip trailing digits
species <- sub("_.*$", "", species)                     # drop anything after underscores if present
species[tip_labels == "reference"] <- "reference"

# last digit from the sample name; blank for 'reference'
last_digit <- ifelse(tip_labels == "reference", "ref",
                     sub("^.*?(\\d)$", "\\1", tip_labels))

tip_df <- data.frame(label = tip_labels,
                     species = species,
                     digit = last_digit,
                     stringsAsFactors = FALSE)

# ---- Read fixed colors from TSV ----
cc <- read_tsv("/xdisk/mcnew/dannyjackson/sulidae/referencelists/colorcodes.tsv", col_types = "cc", progress = FALSE)
# make a named vector like c("BFBO"="#4FE98C", ...)
fixed_map <- setNames(cc$Code, toupper(cc$Species))

# ---- Build species list present in the plots ----
species_levels <- sort(unique(tip_df$species))  # tip_df$species should be like BFBO/PEBO/etc.

# ---- Start palette with NA, fill from file, then handle special/missing ----
pal <- setNames(rep(NA_character_, length(species_levels)), species_levels)

# fill known species from the file
overlap <- intersect(names(fixed_map), species_levels)
pal[overlap] <- fixed_map[overlap]

# always set reference to grey (if present)
if ("reference" %in% names(pal)) pal["reference"] <- "grey55"



# 1) Compute bootstrap support (%) for each internal node on the consensus
bs <- prop.clades(cons_rooted, trees)
bs <- bs[!is.na(bs)]  # remove root NA
bs <- bs / 31 * 100   # rescale if desired

# 2) Define a white→black gradient for values ≥70, and pure white for <70
# bs in percent (0–100), may include NA for the root
node_cols <- rep("white", length(bs))         # default: <70 -> white

ok <- which(!is.na(bs) & bs >= 70)
g  <- (100 - bs[ok]) / 30                     # 70 -> 1 (white), 100 -> 0 (black)
g  <- pmin(pmax(g, 0), 1)                     # clamp to [0,1] just in case
node_cols[ok] <- gray(g)

# Explanation:
# gray(x) expects 0 (black) to 1 (white), so (100-bs)/30 scales 70→100 into 1→0.
# --- Tip label relabeling and coloring by species ---

# 1) Build per-tip colors from the species->color palette
tip_colors <- pal[tip_df$species]
tip_colors[is.na(tip_colors)] <- "black"     # fallback if any species missing in palette

# 2) Relabel: last digit for samples, "ref" for reference
tip_labels_new <- ifelse(tip_df$label == "reference", "ref", tip_df$digit)

# Apply new labels to the plotted tree (safe to modify a copy)
cons_for_plot <- cons_rooted
cons_for_plot$tip.label <- tip_labels_new

# --- Plot with tip text colored by species, and tip symbols in the same color ---
pdf("allbestTree.MRC.pdf", width = 4, height = 10, useDingbats = FALSE)
plot(cons_for_plot,
     type = "phylogram",
     cex = 1.5,
     no.margin = TRUE,
     tip.color = tip_colors,
     label.offset = 0.5)                  # color the numbers

# Add a filled circle at each tip in the species color, black outline
tiplabels(pch = 21, bg = tip_colors, col = "black", cex = 0.9, lwd = 0.6)

# --- Internal node support symbols (from your gradient calc) ---
int_nodes <- (Ntip(cons_for_plot) + 1):(Ntip(cons_for_plot) + cons_for_plot$Nnode)
nodelabels(node = int_nodes, pch = 23, bg = node_cols, col = "black", cex = 1.2, lwd = 0.6)

# Legend for support gradient (unchanged)
legend_vals <- seq(70, 100, by = 5)
legend_cols <- gray((100 - legend_vals) / 30)
legend("bottomleft",
       legend = c("<70", legend_vals),
       pt.bg  = c("white", legend_cols),
       pch    = 21,
       col    = "black",
       pt.cex = 1.0,
       bty    = "n",
       title  = "% Support")

dev.off()



# visualize trees overlapping

# --- Read, root, ladderize (your code) ---
trees <- read.tree("allbestTree.autosomes.nwk")
trees <- lapply(trees, function(x) root(x, outgroup = "reference", resolve.root = TRUE))
trees <- lapply(trees, ladderize)

pdf(file = "densitree.autosomes.pdf", width = 8, height = 10, useDingbats = FALSE)

ggdensitree(trees,
    alpha=.3, colour='steelblue') + 
    geom_tiplab(size=3) + hexpand(.35)

dev.off()

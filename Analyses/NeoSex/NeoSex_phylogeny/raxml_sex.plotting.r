# Plot bestTree

source ~/.bashrc

micromamba create -f /xdisk/mcnew/dannyjackson/sulidae/analyses/raxml_wg/r_trees.yml

micromamba activate r_trees


library(ape)
library(ggplot2)
library(ggtree)
library(treeio)
library(patchwork)
library(scales)   # for hue_pal()
library(readr)
library(stringr)
library(phangorn)  # still fine to keep for other utilities


# Read all trees
tree <- read.tree("trees.CM062600.1.raxml.bestTree")  # or list.files() and lapply(read.tree, ...)
tree_r <- root(tree, outgroup = "reference", resolve.root = TRUE)

# ---- Build tip annotation (shared across trees if same taxa set) ----
tip_labels <- tree_r$tip.label

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


library(treeio)
library(ggtree)
# library(ggplot2)  # already loaded in your session

# Build in clear steps
p <- ggtree(tree_r, size = 0.25) %<+% tip_df +
  geom_tree(linewidth = 0.25) +                           # <-- move size here
  geom_tiplab(
    aes(color = species, label = digit),
    linesize = 0.25,
    geom = "text",
    align = FALSE,
    size = 2,
    linetype = "dotted",
    show.legend = FALSE,
    offset = 0.005
  )  +
  geom_tippoint(aes(color = species), size = 0.5, alpha = 0.9,
                position = position_nudge(x = 0.0015))

p <- p +
  labs(title = "Chromosome 29") +
  scale_color_manual(values = pal, breaks = species_levels, name = "Species") +
  theme(
    plot.title    = element_text(size = 8, hjust = 0.5),
    plot.subtitle = element_text(size = 6, hjust = 0.5, margin = margin(b = 2)),
    plot.margin   = margin(1, 4, 4, 1)
  )

ggsave("tree.png",
       p, width = 3, height = 5, dpi = 300) 

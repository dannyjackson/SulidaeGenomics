suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

infile     <- "bam.abbababa2"  # your ANGSD output
min_sites  <- 50                          # filter noisy blocks (tune)
min_counts <- 3                           # min(ABBA+BABA) per block (tune)

# Pattern sets (match column names without the leading 'X')
ABBA <- c("0110","0220","0330","1001","1221","1331","2002","2112","2332","3003","3113","3223")
BABA <- c("0101","0202","0303","1010","1212","1313","2020","2121","2323","3030","3131","3232")

dt <- fread(infile, sep = "\t", header = TRUE)

# ANGSD headers may be read as X0110, X0101, ...  strip any leading 'X'
stripx <- function(v) sub("^X", "", v)
setnames(dt, names(dt), stripx(names(dt)))

# quick sanity
stopifnot(all(c("CHR","BLOCKstart","BLOCKend","numSites") %in% names(dt)))

# Safe summation helper
sum_cols <- function(DT, cols) {
  cols <- intersect(cols, names(DT))
  if (!length(cols)) return(rep(0, nrow(DT)))
  rowSums(DT[, ..cols], na.rm = TRUE)
}

# Per-block totals and D
dt[, `:=`(
  ABBA_sum  = sum_cols(.SD, ABBA),
  BABA_sum  = sum_cols(.SD, BABA)
)]
dt[, total := ABBA_sum + BABA_sum]
dt[, D     := fifelse(total > 0, (ABBA_sum - BABA_sum) / total, NA_real_)]
dt[, mid   := floor((BLOCKstart + BLOCKend)/2)]

# basic QC filters (optional but recommended)
dtf <- dt[numSites >= min_sites & total >= min_counts & is.finite(D)]

# order chromosomes sensibly if they’re like CM062567.1 …
chr_order <- sort(unique(dtf$CHR), method = "radix")
dtf[, CHR := factor(CHR, levels = chr_order)]

# Plot: one facet per chromosome
p <- ggplot(dtf, aes(x = mid, y = D)) +
  geom_hline(yintercept = 0, linewidth = 0.3, linetype = 2, color = "grey40") +
  geom_point(alpha = 0.5, size = 0.6) +
  # loess smoothing per chr helps see broad trends
  geom_smooth(method = "loess", span = 0.2, se = FALSE, linewidth = 0.8) +
  facet_wrap(~ CHR, scales = "free_x") +
  labs(x = "Position (bp)", y = "D (ABBA−BABA / ABBA+BABA)",
       title = "Block-wise D-statistic from ANGSD doAbbababa2") +
  theme_bw(base_size = 11)

ggsave("D_by_block.faceted.pdf", p, width = 12, height = 8)

# Plot as a manhattan plot 

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

##------------------------------------------------------------
## INPUTS (edit these three as needed)
##------------------------------------------------------------
convfile       <- "/xdisk/mcnew/dannyjackson/sulidae/referencelists/chromconversion.txt"   # two columns mapping: old -> new
# Acceptable column name pairs (case-insensitive):
#   CHR/chrom, scaffold/chromo, old/new, original/chromosome
metric         <- "D"                           # what to plot on y-axis (from your dtf)
metric_cutoff  <- 0                             # horizontal reference line
outdir         <- "out"                         # where to write figure
pop_name       <- "POP"                         # label component for filename
win            <- "blocks"                      # window label component for filename

dir.create(file.path(outdir, "analyses", metric, pop_name, win), recursive = TRUE, showWarnings = FALSE)

##------------------------------------------------------------
## 1) Read/standardize chromosome conversion
##------------------------------------------------------------
conv <- fread(convfile)

colnames(conv) <- c("chromo", "CHR")
# normalize column names to lower, then map to CHR -> chromo
setnames(conv, tolower(names(conv)))
name_map <- list(
  CHR    = intersect(names(conv), c("chr","scaffold","old","original","contig")),
  chromo = intersect(names(conv), c("chromo","chrom","chromosome","new"))
)
# stopifnot(length(name_map$CHR) == 1L, length(name_map$chromo) == 1L)
setnames(conv, name_map$CHR, "CHR")
setnames(conv, name_map$chromo, "chromo")

# coerce target labels to character (e.g., "1","1A","2",...,"Z")
conv[, chromo := as.character(chromo)]

##------------------------------------------------------------
## 2) Start from your filtered ANGSD blocks (dtf) and add mapping
##     Assumes you already produced `dtf` with columns: CHR, mid, D
##------------------------------------------------------------
if (!exists("dtf")) stop("Expected `dtf` from your earlier code (the filtered ABBA/BABA blocks).")

data <- merge(dtf[, .(CHR, position = mid, D)], conv[, .(CHR, chromo)], by = "CHR", all.x = TRUE)

if (anyNA(data$chromo)) {
  missing_scfs <- unique(data[is.na(chromo), CHR])
  stop(sprintf("These CHR IDs were not found in the conversion file:\n%s",
               paste(missing_scfs, collapse = ", ")))
}

##------------------------------------------------------------
## 3) Force chromosome order: 1, 1A, 2:4, 4A, 5:29, Z (drop any not present)
##------------------------------------------------------------
desired_levels <- c("1", "1A", as.character(2:4), "4A", as.character(5:29), "Z")
present_levels <- intersect(desired_levels, unique(data$chromo))
data[, chromo := factor(chromo, levels = present_levels)]

##------------------------------------------------------------
## 4) Build cumulative genome axis (Manhattan-style)
##------------------------------------------------------------
# per-chromosome max coordinate
chr_len <- data[, .(chr_len = max(position, na.rm = TRUE)), by = chromo][order(chromo)]
# cumulative starts (tot) by ordered chromo
chr_len[, tot := shift(cumsum(chr_len), fill = 0)]
# join back and compute BPcum
plot_data <- chr_len[data, on = "chromo"]
plot_data[, BPcum := position + tot]

# x-axis tick centers
axisdf <- plot_data[, .(center = mean(BPcum, na.rm = TRUE)), by = chromo][order(chromo)]

##------------------------------------------------------------
## 5) Colors (alternating two; safe fallback if color1/color2 not defined)
##------------------------------------------------------------
if (!exists("color1")) color1 <- "#bdbdbd"
if (!exists("color2")) color2 <- "#636363"
chrs <- levels(plot_data$chromo)
col_vec <- rep(c(color1, color2), length.out = length(chrs))
names(col_vec) <- chrs

##------------------------------------------------------------
## 6) Plot
##------------------------------------------------------------
yvar <- metric
stopifnot(yvar %in% names(plot_data))  # make sure metric exists

p <- ggplot(plot_data, aes(x = BPcum, y = .data[[yvar]])) +
  geom_point(aes(color = chromo), alpha = 0.8, size = 0.8) +
  scale_color_manual(values = col_vec, guide = "none") +
  scale_x_continuous(breaks = axisdf$center, labels = as.character(axisdf$chromo),
                     guide = guide_axis(n.dodge = 2)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_hline(yintercept = metric_cutoff, linewidth = 0.4) +
  labs(x = "Chromosome", y = yvar,
       title = "Block-wise D across chromosomes") +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

outfile <- file.path(outdir, "analyses", metric, pop_name, win,
                     sprintf("%s.%s.%s.sigline.png", pop_name, metric, win))
ggsave(outfile, p, width = 20, height = 5, units = "in", dpi = 300)

cat("Wrote: ", outfile, "\n")

library(ggplot2)
library(dplyr)
library(stringr)
library(readr)

# Parameters
mu <- 4.6e-09
gen <- 8

# Directory of input files
indir <- "/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/files/autosomes/"
files <- list.files(indir, pattern="*.final.txt", full.names=TRUE)

# Extract data
all_data <- lapply(files, function(f) {
  dat <- read.csv(f, sep="\t")
  
  species <- str_extract(basename(f), "^[A-Z]+")  # grabs BFBO, PEBO, etc.
  indiv <- str_extract(basename(f), "^[A-Z]+[0-9]+")

  time <- (dat$left_time_boundary/mu*gen)
  pop.size <- (1/dat$lambda)/(2*mu)
  
  data.frame(time=time, pop.size=pop.size, species=species, individual=indiv)
})

all_data <- bind_rows(all_data)


# Load color codes
colorcodes <- read_tsv("/xdisk/mcnew/dannyjackson/sulidae/referencelists/colorcodes.tsv")

# Make a named vector for scale_color_manual()
species_colors <- setNames(colorcodes$Code, colorcodes$Species)

# Plot
current_date <- format(Sys.Date(), "%m%d%y")
outfile <- paste0("/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/plots/MSMC_all_species_", current_date, ".autosomes.mingen.pdf")

pdf(outfile, width = 16, height = 6)  # higher resolution and larger canvas

ggplot(all_data, aes(x = time, y = pop.size, color = species, group = individual)) +
  geom_step(linewidth = 1.4, alpha = 0.75) +  # thicker lines
  scale_x_log10(limits = c(1e4, 1e7)) +
  ylim(0, 500000) +
  scale_color_manual(values = species_colors) +  # apply your color scheme
  labs(
    x = "Years before present (log scale)",
    y = "Effective Population Size",
    title = "MSMC across individuals and species"
  ) +
  theme_minimal(base_size = 22) +  # larger base font size
  theme(
    plot.title = element_text(size = 26, hjust = 0.5),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.6),
    legend.position = "null"
  )

dev.off()



# Make a facet plot with one panel per species
suppressMessages({
  library(ggplot2)
  library(scales)
  library(dplyr)
})

# Output dir + date stamp
outdir  <- "/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/plots"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
stamp <- format(Sys.Date(), "%m%d%y")

# Unique species
species_list <- sort(unique(all_data$species))

# Save one PNG per species
for (sp in species_list) {
  df <- filter(all_data, species == sp)

  if (nrow(df) == 0) next

  p_sp <- ggplot(df, aes(x = time, y = pop.size, color = individual, group = individual)) +
  geom_step(linewidth = 1.4, alpha = 0.75) +  # thicker lines
  scale_x_log10(limits = c(1e4, 1e7)) +
  ylim(0, 200000) +
  labs(
    x = "Years before present (log scale)",
    y = "Effective Population Size",
    title = paste0("MSMC per individual — ", sp),
    color = "Individual"
  ) +
  theme_minimal(base_size = 22) +  # larger base font size
  theme(
    plot.title = element_text(size = 26, hjust = 0.5),
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.6),
  )

  outfile <- file.path(outdir, paste0("MSMC_", sp, "_", stamp, ".mingen.pdf"))
  ggsave(outfile, p_sp, width = 10, height = 7, dpi = 300)
  message("Saved: ", outfile)
}
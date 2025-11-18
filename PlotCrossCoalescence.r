# BRBO_COBO
# PEBO_NABO
# PEBO_MABO
# BFBO_NABO
# BFBO_MABO
# BFBO_PEBO
# MABO_NABO # results/MABO305_NABO402_.msmc2.final.txt

# Parameters
mu  <- 1.913e-9
gen <- 14.11

# Output file
pdf("MSMC2_cross_coalescence.pdf", width = 7, height = 5)

# List all result files
files <- list.files("results", pattern = "\\.msmc2\\.final\\.txt$", full.names = TRUE)

# Set up empty plot
plot(NA, NA,
     xlim = c(1000000, 10000000),
     ylim = c(0, 1),
     xlab = "Years ago",
     ylab = "Relative cross-coalescence rate (RCCR)",
     main = "MSMC2 cross-coalescence curves")

# Loop through and plot all in same color
for (f in files) {
  dat <- read.table(f, header = TRUE)
  years <- dat$left_time_boundary / mu * gen
  rccr  <- 2 * dat$lambda_01 / (dat$lambda_00 + dat$lambda_11)
  
  lines(years, rccr, type = "s", col = "#4FE98C", lwd = 1)
  # lines(years, rccr, type = "s", col = "#B32904", lwd = 1)
}


dev.off()


####################################################
# Color by first reference taxa
# Parameters
mu  <- 1.913e-9
gen <- 14.11

# Output file
pdf("MSMC2_cross_coalescence_byBFBO.pdf", width = 7, height = 5)

# List all result files
files <- list.files("results", pattern = "\\.msmc2\\.final\\.txt$", full.names = TRUE)

# Extract BFBO ID from each filename (e.g., BFBO501)
bfbo_ids <- sub(".*(BFBO[0-9]+).*", "\\1", basename(files))

# Assign colors by BFBO individual
unique_bfbo <- unique(bfbo_ids)
bfbo_cols <- setNames(rainbow(length(unique_bfbo)), unique_bfbo)

# Set up empty plot
plot(NA, NA,
     xlim = c(1e6, 1e7),
     ylim = c(0, 1),
     xlab = "Years ago",
     ylab = "Relative cross-coalescence rate (RCCR)",
     main = "MSMC2 cross-coalescence curves by BFBO individual")

# Loop through files and plot colored by BFBO
for (f in files) {
  dat <- read.table(f, header = TRUE)
  years <- dat$left_time_boundary / mu * gen
  rccr  <- 2 * dat$lambda_01 / (dat$lambda_00 + dat$lambda_11)
  
  bfbo_id <- sub(".*(BFBO[0-9]+).*", "\\1", basename(f))
  col_use <- bfbo_cols[bfbo_id]
  
  lines(years, rccr, type = "s", col = col_use, lwd = 1.5)
}

# Add legend
legend("topright",
       legend = unique_bfbo,
       col = bfbo_cols,
       lwd = 2, cex = 0.8, bg = "white",
       title = "BFBO individual")

dev.off()

####################################################
# Color by second reference taxa
# Parameters
# Parameters
mu  <- 1.913e-9
gen <- 14.11

# Output file
pdf("MSMC2_cross_coalescence_byBRBO.pdf", width = 7, height = 5)

# List all result files
files <- list.files("results", pattern = "\\.msmc2\\.final\\.txt$", full.names = TRUE)

# Extract PEBO ID from each filename (e.g., PEBO601)
pebo_ids <- sub(".*?(BRBO[0-9]+).*", "\\1", basename(files))

# Assign colors by PEBO individual
unique_pebo <- unique(pebo_ids)
pebo_cols <- setNames(rainbow(length(unique_pebo)), unique_pebo)

# Set up empty plot
plot(NA, NA,
     xlim = c(1e6, 1e7),
     ylim = c(0, 1),
     xlab = "Years ago",
     ylab = "Relative cross-coalescence rate (RCCR)",
     main = "MSMC2 cross-coalescence curves by PEBO individual")

# Loop through files and plot colored by PEBO
for (f in files) {
  dat <- read.table(f, header = TRUE)
  years <- dat$left_time_boundary / mu * gen
  rccr  <- 2 * dat$lambda_01 / (dat$lambda_00 + dat$lambda_11)
  
  pebo_id <- sub(".*?(BRBO[0-9]+).*", "\\1", basename(files))
  col_use <- pebo_cols[pebo_id]
  
  lines(years, rccr, type = "s", col = col_use, lwd = 1.5)
}

# Add legend
legend("topright",
       legend = unique_pebo,
       col = pebo_cols,
       lwd = 2, cex = 0.8, bg = "white",
       title = "BRBO individual")

dev.off()

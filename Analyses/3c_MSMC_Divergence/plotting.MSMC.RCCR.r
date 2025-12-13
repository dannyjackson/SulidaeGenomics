####################################################
# General MSMC2 RCCR plotting by "side" of pair
# - Filenames like: SPECIES1ID_SPECIES2ID.msmc2.final.txt
#   e.g., BFBO501_PEBO601.msmc2.final.txt
# - Color by first individual (left) or second (right)
####################################################

# Params you may tweak
mu  <- 1.913e-9
gen <- 14.11
results_dir <- "results"  # folder containing *.msmc2.final.txt
xlim_years  <- c(1e6, 1e7) # or set to NULL for auto
ylim_rccr   <- c(0, 1)     # RCCR is [0,1]

# Helper to make distinct colors
distinct_cols <- function(n) grDevices::rainbow(n, v = 0.9, s = 0.9)

# Parse a filename into left/right IDs (e.g., "BFBO501", "PEBO601")
parse_ids <- function(fn_base) {
  # Expect something like ABCD1234_WXYZ567.msmc2.final.txt
  m <- regexec("^([A-Z]{4}\\d+)_([A-Z]{4}\\d+)\\.msmc2\\.final\\.txt$", fn_base)
  hits <- regmatches(fn_base, m)
  if (length(hits) == 0 || length(hits[[1]]) < 3) return(c(NA_character_, NA_character_))
  hits[[1]][2:3]
}

# Plot function
plot_rccr_by_side <- function(side = c("left", "right"),
                              out_pdf = NULL) {
  side <- match.arg(side)

  # Gather files
  files <- list.files(results_dir, pattern = "\\.msmc2\\.final\\.txt$", full.names = TRUE)
  if (length(files) == 0) {
    stop("No *.msmc2.final.txt files found in '", results_dir, "'.")
  }

  # Extract IDs per file
  bases <- basename(files)
  ids   <- lapply(bases, parse_ids)
  ids_mat <- do.call(rbind, ids)
  colnames(ids_mat) <- c("left_id", "right_id")

  # Drop any non-matching files
  keep <- !is.na(ids_mat[,1]) & !is.na(ids_mat[,2])
  files <- files[keep]
  bases <- bases[keep]
  ids_mat <- ids_mat[keep, , drop = FALSE]
  if (length(files) == 0) stop("No files matched the expected pattern.")

  # Choose which column to color by
  color_ids <- if (side == "left") ids_mat[, "left_id"] else ids_mat[, "right_id"]
  unique_ids <- unique(color_ids)
  cols <- setNames(distinct_cols(length(unique_ids)), unique_ids)

  # Title / output
  # Use species codes (first 4 letters) from the *first* file to build a title
  # This is just cosmetic; we’ll still parse every file individually
  sp_left  <- sub("\\d+", "", ids_mat[1, "left_id"])
  sp_right <- sub("\\d+", "", ids_mat[1, "right_id"])
  who <- if (side == "left") paste0(sp_left, " individual") else paste0(sp_right, " individual")
  main_title <- sprintf("MSMC2 cross-coalescence colored by %s", who)

  # Output pdf name if not provided
  if (is.null(out_pdf)) {
    out_pdf <- sprintf("MSMC2_RCCR_by_%s.pdf",
                       if (side == "left") sp_left else sp_right)
  }

  # Read one file to get axis ranges if auto
  # (Alternatively, compute over all files; this is usually fine)
  compute_ranges <- function(ff) {
    dat <- read.table(ff, header = TRUE)
    years <- dat$left_time_boundary / mu * gen
    rccr  <- 2 * dat$lambda_01 / (dat$lambda_00 + dat$lambda_11)
    list(x = range(years, finite = TRUE), y = range(rccr, finite = TRUE))
  }
  if (is.null(xlim_years) || is.null(ylim_rccr)) {
    xr <- c(Inf, -Inf)
    yr <- c(Inf, -Inf)
    for (f in files) {
      rr <- compute_ranges(f)
      xr <- c(min(xr[1], rr$x[1], na.rm = TRUE), max(xr[2], rr$x[2], na.rm = TRUE))
      yr <- c(min(yr[1], rr$y[1], na.rm = TRUE), max(yr[2], rr$y[2], na.rm = TRUE))
    }
    if (is.null(xlim_years)) xlim_years <- xr
    if (is.null(ylim_rccr))  ylim_rccr  <- yr
  }

  # Plot
  pdf(out_pdf, width = 7, height = 5)
  on.exit(dev.off(), add = TRUE)

  plot(NA, NA,
       xlim = xlim_years,
       ylim = ylim_rccr,
       xlab = "Years ago",
       ylab = "Relative cross-coalescence rate (RCCR)",
       main = main_title)

  # Draw curves
  for (i in seq_along(files)) {
    f <- files[i]
    dat <- read.table(f, header = TRUE)
    years <- dat$left_time_boundary / mu * gen
    rccr  <- 2 * dat$lambda_01 / (dat$lambda_00 + dat$lambda_11)

    id_use <- color_ids[i]
    col_use <- cols[[id_use]]

    # Step-like lines as in MSMC plots
    lines(years, rccr, type = "s", col = col_use, lwd = 1.5)
  }

  # Legend
  legend("topright",
         legend = names(cols),
         col = unname(cols),
         lwd = 2, cex = 0.8, bg = "white",
         title = if (side == "left") paste0(sp_left, " individual")
                 else paste0(sp_right, " individual"))
  message("Wrote: ", out_pdf)
}

# ---- Run both variants (by first/left and by second/right) ----
plot_rccr_by_side(side = "left")   # colors by first individual in each pair
plot_rccr_by_side(side = "right")  # colors by second individual in each pair

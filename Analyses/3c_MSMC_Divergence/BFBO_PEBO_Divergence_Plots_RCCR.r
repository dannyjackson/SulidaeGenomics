library(ggplot2)
library(dplyr)
library(readr)

setwd("/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/files/species/inferring_split_times")

data <- read.table("BFBO_PEBO.final.txt", header = TRUE)

mu <- 1.913e-9

# Compute RCCR once (independent of generation time)
df_base <- data %>%
  mutate(
    CCR = 2 * lambda_01 / (lambda_00 + lambda_11),
    CCR = pmin(pmax(CCR, 0), 1)   # clamp to [0,1] for display
  )

# Make three versions with different generation times (8, 13, and 18 years)
df_8 <- df_base %>%
  mutate(
    gen = 8,
    time_years = (left_time_boundary / mu) * gen
  )

df_13 <- df_base %>%
  mutate(
    gen = 13,
    time_years = (left_time_boundary / mu) * gen
  )

df_18 <- df_base %>%
  mutate(
    gen = 18,
    time_years = (left_time_boundary / mu) * gen
  )

df_plot <- bind_rows(df_8, df_13, df_18) %>%
  mutate(
    gen = factor(gen, levels = c(8, 13, 18))
  )


# ----------------------------
# Temperature quantiles
# temperatures.csv columns like:
# Time (kyr BP), 2.5%, 5%, 25%, 50%, 75%, 95%, 97.5%
# ----------------------------
temps <- read_csv("temperatures.csv", show_col_types = FALSE)


# Find the time column robustly and convert to years BP
time_col <- names(temps)[grepl("^Time", names(temps), ignore.case = TRUE)][1]

temps <- temps %>%
  rename(Time_kyr_BP = all_of(time_col)) %>%
  mutate(
    time_years = Time_kyr_BP * 1000
  )


# Map quantile columns (after renaming: 2_5pct, 50pct, 97_5pct, etc.)
q2_5  <- names(temps)[grepl("2.5",  names(temps))]
q50   <- names(temps)[grepl("^50$", names(temps))]
q97_5 <- names(temps)[grepl("97.5", names(temps))]

stopifnot(length(q2_5) == 1, length(q50) == 1, length(q97_5) == 1)

temps <- temps %>%
  transmute(
    time_years,
    t_lo = .data[[q2_5]],
    t_md = .data[[q50]],
    t_hi = .data[[q97_5]]
  )


# ----------------------------
# Overlay plot with secondary y-axis
# (RCCR = left axis; ΔGAST = right axis)
# ----------------------------
# We need a linear transform to put temperature on the RCCR (0..1) scale.
# Use the 2.5–97.5% range within the plotting window to define the mapping.
xmax <- 1e6

temps_win <- temps %>% filter(time_years >= 0, time_years <= xmax)

t_min <- min(temps_win$t_lo, na.rm = TRUE)
t_max <- max(temps_win$t_hi, na.rm = TRUE)

# map: temp -> y_on_rccr = (temp - t_min)/(t_max - t_min)
# inverse for sec.axis: y -> temp = y*(t_max-t_min) + t_min
p_all <- ggplot(
  df_plot,
  aes(
    x = time_years,
    y = CCR,
    linetype = gen,
    linewidth = gen
  )
) +
  # temperature ribbon + median (mapped onto 0..1)
  geom_ribbon(
    data = temps_win,
    aes(
      x = time_years,
      ymin = (t_lo - t_min) / (t_max - t_min),
      ymax = (t_hi - t_min) / (t_max - t_min)
    ),
    inherit.aes = FALSE,
    alpha = 0.25
  ) +
  geom_line(
    data = temps_win,
    aes(
      x = time_years,
      y = (t_md - t_min) / (t_max - t_min)
    ),
    inherit.aes = FALSE,
    linewidth = 0.2
  ) +
  # RCCR curves
  geom_step(color = "#7570b3") +
  scale_linetype_manual(values = c(`8` = "dashed", `13` = "solid", `18` = "dotted")) +
  scale_linewidth_manual(values = c(`8` = 0.5, `13` = 0.9, `18` = 0.5)) +
  coord_cartesian(xlim = c(0, xmax), ylim = c(0, 1)) +
  scale_y_continuous(
    name = "Relative cross-coalescence rate (RCCR)",
    limits = c(0, 1),
    sec.axis = sec_axis(
      trans = ~ . * (t_max - t_min) + t_min,
      name = "Δ Global average surface temperature (°C; relative to 0–5 ka)"
    )
  ) +
  labs(
    x = "Years before present",
    linetype = "Generation time (years)",
    linewidth = "Generation time (years)"
  ) +
  theme_classic()

pdf("BFBO_PEBO_RCCR_with_GAST_quantiles.pdf", width = 14, height = 5)
print(p_all)
dev.off()


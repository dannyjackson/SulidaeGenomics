library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)

setwd("/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/files/species/inferring_split_times")

mu <- 1.91e-9
xmax <- 50000

LGM_time <- 21000
LGM_color <- "#9d196aff"   # soft pink

# ----------------------------
# Helper to build RCCR df for a pair
# ----------------------------
make_rccr_df <- function(file, mu) {
  data <- read.table(file, header = TRUE)

  df_base <- data %>%
    mutate(
      CCR = 2 * lambda_01 / (lambda_00 + lambda_11),
      CCR = pmin(pmax(CCR, 0), 1)
    )

  bind_rows(
    df_base %>% mutate(gen = 8,  time_years = (left_time_boundary / mu) * gen),
    df_base %>% mutate(gen = 13, time_years = (left_time_boundary / mu) * gen),
    df_base %>% mutate(gen = 18, time_years = (left_time_boundary / mu) * gen)
  ) %>%
    mutate(gen = factor(gen, levels = c(8, 13, 18)))
}

# ----------------------------
# Temperature quantiles
# ----------------------------
temps <- read_csv("temperatures.csv", show_col_types = FALSE)

time_col <- names(temps)[grepl("^Time", names(temps), ignore.case = TRUE)][1]
temps <- temps %>%
  rename(Time_kyr_BP = all_of(time_col)) %>%
  mutate(time_years = Time_kyr_BP * 1000)

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

temps_win <- temps %>% filter(time_years >= 0, time_years <= xmax)
t_min <- min(temps_win$t_lo, na.rm = TRUE)
t_max <- max(temps_win$t_hi, na.rm = TRUE)

# ----------------------------
# Panel 1: temperature only
# ----------------------------
# Choose nice temperature breaks (edit if you want fewer/more)
temp_breaks <- pretty(c(t_min, t_max), n = 5)

# Convert temp breaks → RCCR scale (0–1)
temp_breaks_scaled <- (temp_breaks - t_min) / (t_max - t_min)

p_temp <- ggplot(temps_win, aes(x = time_years)) +
  geom_ribbon(
    aes(
      ymin = (t_lo - t_min) / (t_max - t_min),
      ymax = (t_hi - t_min) / (t_max - t_min)
    ),
    alpha = 0.25
  ) +
  geom_line(
    aes(y = (t_md - t_min) / (t_max - t_min)),
    linewidth = 0.25
  ) +
  coord_cartesian(xlim = c(0, xmax), ylim = c(0, 1)) +
  scale_y_continuous(
    name = "Temp (°C)",
    breaks = temp_breaks_scaled,
    labels = round(temp_breaks, 1)
  ) +
  labs(x = NULL,
    title = "Δ Global avg. surface temperature (relative to 0–5 ka)") +
  geom_vline(
  xintercept = LGM_time,
  linewidth = 2.5,
  color = LGM_color,
  alpha = 0.35
  ) +
  theme_classic()


# ----------------------------
# Panel 2: BFBO–PEBO (your original styling)
# ----------------------------
df_bfbo_pebo <- make_rccr_df("BFBO_PEBO.final.txt", mu)

p_bfbo_pebo <- ggplot(
  df_bfbo_pebo,
  aes(x = time_years, y = CCR, linetype = gen, linewidth = gen)
) +
  geom_step(color = "#7570b3") +
  scale_linetype_manual(values = c(`8` = "dashed", `13` = "solid", `18` = "dotted")) +
  scale_linewidth_manual(values = c(`8` = 0.5, `13` = 0.9, `18` = 0.5)) +
  coord_cartesian(xlim = c(0, xmax), ylim = c(0, 1)) +
  labs(
    x = NULL,
    y = "RCCR",
    linetype = "Generation time (years)",
    linewidth = "Generation time (years)",
    title = "Blue-footed vs. Peruvian"
  ) +
  geom_vline(
  xintercept = LGM_time,
  linewidth = 2.5,
  color = LGM_color,
  alpha = 0.35
  ) +
  theme(
  legend.justification = c(0, 0.5),
  legend.title.align = 0
  ) +
  theme_classic()

# ----------------------------
# Panel 3: MABO–NABO (same as BFBO–PEBO)
# ----------------------------
df_mabo_nabo <- make_rccr_df("MABO_NABO.final.txt", mu)

p_mabo_nabo <- ggplot(
  df_mabo_nabo,
  aes(x = time_years, y = CCR, linetype = gen, linewidth = gen)
) +
  geom_step(color = "#7570b3") +
  scale_linetype_manual(values = c(`8` = "dashed", `13` = "solid", `18` = "dotted")) +
  scale_linewidth_manual(values = c(`8` = 0.5, `13` = 0.9, `18` = 0.5)) +
  coord_cartesian(xlim = c(0, xmax), ylim = c(0, 1)) +
  labs(
    x = "Years before present",
    y = "RCCR",
    title = "Masked vs. Nazca"
  ) +
  geom_vline(
  xintercept = LGM_time,
  linewidth = 2.5,
  color = LGM_color,
  alpha = 0.35
  ) +
  theme_classic()  +
  geom_point(
  aes(x = 5000000, y = 5, color = "LGM"),
  shape = 15,   # filled square
  size = 4,
  alpha = 0.35
) +
  scale_color_manual(
    name = "Glacial events",
    values = c(
      "LGM" = LGM_color)
  ) +
  guides(
    linetype = "none",
    linewidth = "none"
  ) +
  theme(
  legend.justification = c(0, 0.5),
  legend.title.align = 0
  )


# ----------------------------
# Stack into one figure
# ----------------------------
p_stack <- p_temp / p_bfbo_pebo / p_mabo_nabo +
  plot_layout(heights = c(0.8, 1.1, 1.1))

ggsave("RCCR_temp_BFBO_PEBO_MABO_NABO_stacked.50k.png", p_stack, width = 7, height = 5)

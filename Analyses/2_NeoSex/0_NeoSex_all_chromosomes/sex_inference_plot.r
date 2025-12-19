# Load packages
library(tidyverse)

# --- Input files ---
z_file   <- "Z_matrix.csv"
sex_file  <- "/xdisk/mcnew/dannyjackson/sulidae/referencelists/sexID.csv"
conv_file <- "/xdisk/mcnew/dannyjackson/sulidae/referencelists/chromconversion.txt"

# --- Read data ---
z  <- read_tsv(z_file, show_col_types = FALSE)
sex  <- read_csv(sex_file, col_names = c("individual", "sex"), show_col_types = FALSE)
conv <- read_csv(conv_file, col_names = c("chrom_num", "chrom"), show_col_types = FALSE)

colnames(z) = c("chrom", "individual", "H_O", "norm_depth")

z <- z %>%
  left_join(sex, by = "individual") %>%
  mutate(sex = fct_explicit_na(sex, na_level = "Unknown"))

# --- Long format for 4 panels in requested order ---
plot_df <- z %>%
  select(individual, sex, chrom, H_O, norm_depth) %>%
  pivot_longer(cols = c(H_O, norm_depth), names_to = "metric", values_to = "value") %>%
  mutate(
    panel = case_when(
      chrom == "Z1" & metric == "H_O"        ~ "Z1 Het",
      chrom == "Z1" & metric == "norm_depth" ~ "Z1 Depth",
      chrom == "Z2" & metric == "H_O"        ~ "Z2 Het",
      chrom == "Z2" & metric == "norm_depth" ~ "Z2 Depth"
    ),
    panel = factor(panel, levels = c("Z1 Het","Z2 Het","Z1 Depth","Z2 Depth"))
  ) %>%
  filter(!is.na(panel))

# x positions for the 4 panels so we can draw connecting lines
plot_df <- plot_df %>%
  mutate(x = as.integer(panel))

# --- Colors ---
sex_cols <- c(Female="#E26D5A", Male="#4F7CAC", Unknown="grey70")

# --- Plot: points + individual lines across panels ---
gg <- ggplot(plot_df, aes(x = x, y = value)) +
  # connect each individual across the 4 panels
  geom_line(aes(group = individual, color = sex), alpha = 0.35, linewidth = 0.35) +
  geom_point(aes(color = sex), alpha = 0.75, size = 1.2) +
  scale_color_manual(values = sex_cols) +
  scale_x_continuous(
    breaks = 1:4,
    labels = levels(plot_df$panel),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 18) +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(size = 13),
    panel.grid.minor = element_blank()
  )

# --- Save ---
ggsave("sex_inference.pdf", gg, width = 5, height = 4, dpi = 300)

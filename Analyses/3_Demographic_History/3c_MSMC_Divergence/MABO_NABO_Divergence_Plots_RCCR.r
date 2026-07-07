library(ggplot2)
library(dplyr)

setwd("/xdisk/mcnew/dannyjackson/sulidae/analyses/msmc/files/species/inferring_split_times")

data <- read.table("MABO_NABO.final.txt", header = TRUE)

mu <- 1.91e-9

# Compute RCCR once (independent of generation time)
df_base <- data %>%
  mutate(
    CCR = 2 * lambda_01 / (lambda_00 + lambda_11),
    CCR = pmin(pmax(CCR, 0), 1)   # clamp to [0,1] for display
  )

# Make two versions with different generation times (8 and 18 years)
df_8 <- df_base %>%
  mutate(
    gen = 8,
    time_years = (left_time_boundary / mu) * gen
  )

df_18 <- df_base %>%
  mutate(
    gen = 18,
    time_years = (left_time_boundary / mu) * gen
  )

df_plot <- bind_rows(df_8, df_18) %>%
  mutate(
    gen = factor(gen, levels = c(8, 18))
  )

# ---------- Plot 0–5,000,000 years ----------
p_all <- ggplot(df_plot, aes(x = time_years, y = CCR, linetype = gen)) +
  geom_step(color = "#7570b3", linewidth = 0.9) +
  scale_linetype_manual(values = c(`8` = "dashed", `18` = "solid")) +
  coord_cartesian(xlim = c(0, 2000000), ylim = c(0, 1)) +
  labs(
    x = "Years before present",
    y = "Relative cross-coalescence rate (RCCR)",
    linetype = "Generation time (years)"
  ) +
  theme_classic()

pdf("MABO_NABO_RCCR_gen8_vs_gen18.pdf", width = 7, height = 5)
print(p_all)
dev.off()

# ---------- Plot 0–1,000,000 years (zoomed) ----------
p_zoom <- ggplot(df_plot, aes(x = time_years, y = CCR, linetype = gen)) +
  geom_step(color = "#7570b3", linewidth = 0.9) +
  scale_linetype_manual(values = c(`8` = "dashed", `18` = "solid")) +
  coord_cartesian(xlim = c(0, 1000000), ylim = c(0, 1)) +
  labs(
    x = "Years before present",
    y = "Relative cross-coalescence rate (RCCR)",
    linetype = "Generation time (years)"
  ) +
  theme_classic()

pdf("MABO_NABO_RCCR_gen8_vs_gen18.zoomed.pdf", width = 7, height = 5)
print(p_zoom)
dev.off()



# ---------- Plot 0–1,000,000 years (zoomed) ----------
p_zoom_zoom <- ggplot(df_plot, aes(x = time_years, y = CCR, linetype = gen)) +
  geom_step(color = "#7570b3", linewidth = 0.9) +
  scale_linetype_manual(values = c(`8` = "dashed", `18` = "solid")) +
  coord_cartesian(xlim = c(0, 250000), ylim = c(0, 1)) +
  labs(
    x = "Years before present",
    y = "Relative cross-coalescence rate (RCCR)",
    linetype = "Generation time (years)"
  ) +
  theme_classic()

pdf("MABO_NABO_RCCR_gen8_vs_gen18.zoomzoomed.pdf", width = 7, height = 5)
print(p_zoom_zoom)
dev.off()





#### Plot the average gen time
library(ggplot2)
library(dplyr)

data <- read.table("MABO_NABO.final.txt", header = TRUE)

mu  <- 1.91e-9
gen <- 13

df <- data %>%
  mutate(
    CCR = 2 * lambda_01 / (lambda_00 + lambda_11),
    CCR = pmin(pmax(CCR, 0), 1),        # clamp for display
    time_years = (left_time_boundary / mu) * gen
  )

p_all <- ggplot(df, aes(x = time_years, y = CCR)) +
  geom_step(color = "#7570b3", linewidth = 0.9) +
  coord_cartesian(xlim = c(0, 2000000), ylim = c(0, 1)) +
  labs(
    x = "Years before present",
    y = "Relative cross-coalescence rate (RCCR)"
  ) +
  theme_classic()

pdf("MABO_NABO_RCCR_gen13.pdf", width = 7, height = 5)
print(p_all)
dev.off()

p_zoom <- ggplot(df, aes(x = time_years, y = CCR)) +
  geom_step(color = "#7570b3", linewidth = 0.9) +
  coord_cartesian(xlim = c(0, 1000000), ylim = c(0, 1)) +
  labs(
    x = "Years before present",
    y = "Relative cross-coalescence rate (RCCR)"
  ) +
  theme_classic()

pdf("MABO_NABO_RCCR_gen13.zoomed.pdf", width = 7, height = 5)
print(p_zoom)
dev.off()

p_zoom_zoom <- ggplot(df, aes(x = time_years, y = CCR)) +
  geom_step(color = "#7570b3", linewidth = 0.9) +
  coord_cartesian(xlim = c(0, 250000), ylim = c(0, 1)) +
  labs(
    x = "Years before present",
    y = "Relative cross-coalescence rate (RCCR)"
  ) +
  theme_classic()

pdf("MABO_NABO_RCCR_gen13.zoomzoomed.pdf", width = 7, height = 5)
print(p_zoom_zoom)
dev.off()



# RCCR already computed as `ccr`
df <- df %>% arrange(time_years)

cross_idx <- which(df$CCR < 0.5)
t_cross <- df$time_years[cross_idx]

t_cross
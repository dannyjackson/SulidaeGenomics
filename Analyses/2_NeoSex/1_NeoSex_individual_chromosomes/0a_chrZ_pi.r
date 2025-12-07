library(readr); library(dplyr); library(ggplot2)

win <- read_tsv("CM062600.1.pi_1000000.windowed.pi",
                show_col_types = FALSE,
                col_types = cols(
                  CHROM = col_character(),
                  BIN_START = col_double(),
                  BIN_END = col_double(),
                  N_VARIANTS = col_double(),
                  PI = col_double()
                ))

ggplot(win, aes(x = (BIN_START + BIN_END)/2, y = PI)) +
  geom_line(linewidth = 0.4) +
  labs(title = "Windowed nucleotide diversity (π)",
       x = "Position (bp)", y = "π per 1 Mb") +
  theme_bw()


ggsave("pi.CM062600.1.allindv.png", width = 12, height = 3, dpi = 300)

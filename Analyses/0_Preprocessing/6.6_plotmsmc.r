#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
OUTDIR = args[1]
IND = args[2]

data <- read.csv(paste0(OUTDIR, "/output/msmc_output.msmc_", IND, "_012525.final.txt"), sep ='\t')

mu <- 2.3e-09
gen <- 7

time <- (data$left_time_boundary/mu*gen)
pop.size <- (1/data$lambda)/(2*mu)

current_date <- format(Sys.Date(), "%m%d%y")

png(paste0(OUTDIR, "/plots/msmc.012525.", IND, ".", current_date, ".cropped.png"), width=500, height=500)

plot(time, pop.size, type="s", 
    xlab="log Years before present", 
    ylab="Effective Population Size", 
    log="x", ylim = c(0,500000))
##########################################

dev.off()


png(paste0(OUTDIR, "/plots/msmc.012525.", IND, ".", current_date, ".png"), width=500, height=500)

plot(time, pop.size, type="s", 
    xlab="log Years before present", 
    ylab="Effective Population Size", 
    log="x")
##########################################

dev.off()


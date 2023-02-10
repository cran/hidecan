## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(tibble)
library(dplyr)
library(purrr)
library(stringr)

## ----setup--------------------------------------------------------------------
library(hidecan)

## ----get-example-data---------------------------------------------------------
x <- get_example_data()

str(x, max.level = 1)

## ----example-data-gwas--------------------------------------------------------
head(x[["GWAS"]])

## ----example-data-de----------------------------------------------------------
head(x[["DE"]])

## ----example-data-can---------------------------------------------------------
head(x[["CAN"]])

## ----hidecan-plot, fig.width = 10, fig.height = 10----------------------------
hidecan_plot(
  gwas_list = x[["GWAS"]],          ## data-frame of GWAS results          
  de_list = x[["DE"]],              ## data-frame of DE results              
  can_list = x[["CAN"]],            ## data-frame of candidate genes            
  score_thr_gwas = -log10(0.0001),  ## sign. threshold for GWAS
  score_thr_de = -log10(0.05),      ## sign. threshold for DE
  log2fc_thr = 0                    ## log2FC threshold for DE
)

## ----hidecan-plot-gwas-can-only, fig.width = 10, fig.height = 10--------------
hidecan_plot(
  gwas_list = x[["GWAS"]],          
  can_list = x[["CAN"]],            
  score_thr_gwas = 4
)

## ----hidecan-plot-with-empty-chrom, fig.width = 10, fig.height = 10-----------
## Chromosomes 0, 6, 9 and 10 are empty
hidecan_plot(
  gwas_list = x[["GWAS"]],          
  can_list = x[["CAN"]],            
  score_thr_gwas = 5
)

## ----hidecan-plot-without-empty-chrom, fig.width = 10, fig.height = 8---------
hidecan_plot(
  gwas_list = x[["GWAS"]],          
  can_list = x[["CAN"]],            
  score_thr_gwas = 5,
  remove_empty_chrom = TRUE
)

## ----hidecan-select-chroms, fig.width = 10, fig.height = 3--------------------
hidecan_plot(
  gwas_list = x[["GWAS"]],                    
  de_list = x[["DE"]],                          
  can_list = x[["CAN"]],                  
  score_thr_gwas = -log10(0.0001),  
  score_thr_de = -log10(0.05),      
  log2fc_thr = 0,
  chroms = c("ST4.03ch07", "ST4.03ch08")
)

## ----hidecan-chrom-limits-all-chrom, fig.width = 10, fig.height = 10----------
hidecan_plot(
  gwas_list = x[["GWAS"]],                    
  de_list = x[["DE"]],                          
  can_list = x[["CAN"]],                  
  score_thr_gwas = -log10(0.0001),  
  score_thr_de = -log10(0.05),      
  log2fc_thr = 0,
  chrom_limits = c(10e6, 20e6)
)

## ----hidecan-chrom-limits-some-chrom, fig.width = 10, fig.height = 10---------
hidecan_plot(
  gwas_list = x[["GWAS"]],                    
  de_list = x[["DE"]],                          
  can_list = x[["CAN"]],                  
  score_thr_gwas = -log10(0.0001),  
  score_thr_de = -log10(0.05),      
  log2fc_thr = 0,
  chrom_limits = list("ST4.03ch01" = c(10e6, 20e6),
                      "ST4.03ch05" = c(30e6, 40e6))
)

## ----hidecan-select-chroms-and-chrom-lims, fig.width = 10, fig.height = 3-----
hidecan_plot(
  gwas_list = x[["GWAS"]],                    
  de_list = x[["DE"]],                          
  can_list = x[["CAN"]],                  
  score_thr_gwas = -log10(0.0001),  
  score_thr_de = -log10(0.05),      
  log2fc_thr = 0,
  chroms = c("ST4.03ch07", "ST4.03ch08"),
  chrom_limits = list("ST4.03ch07" = c(50e6, 55e6),
                      "ST4.03ch08" = c(45e6, 50e6))
)

## ----hidecan-genes-colour-log2fc, fig.width = 10, fig.height = 10-------------
hidecan_plot(
  gwas_list = x[["GWAS"]],          
  de_list = x[["DE"]],              
  can_list = x[["CAN"]],            
  score_thr_gwas = -log10(0.0001),
  score_thr_de = -log10(0.05),
  log2fc_thr = 0,
  colour_genes_by_score = FALSE
)

## ----making-small-dataset-----------------------------------------------------
library(dplyr)
library(purrr)
library(stringr)

## Retaining only markers and genes on chromosomes 7 and 8
x_small <- x |> 
  map(~ filter(.x, str_detect(chromosome, "(07|08)")))

## ----creating-second-gwas-tibble----------------------------------------------
## Creating a second GWAS result tibble by shuffling 
## the marker scores from the original data
gwas_1 <- x_small[["GWAS"]]
gwas_2 <- gwas_1 |> 
  mutate(score = sample(score))

## ----hidecan-multiple-gwas-input, fig.width = 10, fig.height = 3--------------
hidecan_plot(
  gwas_list = list(gwas_1, gwas_2),
  score_thr_gwas = -log10(0.0001),
  score_thr_de = -log10(0.05),
  log2fc_thr = 0
)

## ----multiple-gwas-input-named, fig.width = 10, fig.height = 3----------------
hidecan_plot(
  gwas_list = list("Trait 1" = gwas_1, 
                   "Trait 2" = gwas_2),
  score_thr_gwas = -log10(0.0001),
  score_thr_de = -log10(0.05),
  log2fc_thr = 0
)

## ----make-potato-chrom-length-------------------------------------------------
library(tibble)

## Chromosomes length as recorded in Ensembl Plants
potato_chrom_length <- c(
  ST4.03ch00 = 45813526,
  ST4.03ch01 = 88663952,
  ST4.03ch02 = 48614681,
  ST4.03ch03 = 62190286,
  ST4.03ch04 = 72208621,
  ST4.03ch05 = 52070158,
  ST4.03ch06 = 59532096,
  ST4.03ch07 = 56760843,
  ST4.03ch08 = 56938457,
  ST4.03ch09 = 61540751,
  ST4.03ch10 = 59756223,
  ST4.03ch11 = 45475667,
  ST4.03ch12 = 61165649
) |> 
  ## turn a named vector into a tibble
  enframe(name = "chromosome",
          value = "length")

head(potato_chrom_length)

## ----hidecan-chrom-length, fig.width = 10, fig.height = 10--------------------
hidecan_plot(
  gwas_list = x[["GWAS"]],          
  de_list = x[["DE"]],              
  can_list = x[["CAN"]],
  score_thr_gwas = -log10(0.0001),
  score_thr_de = -log10(0.05),
  log2fc_thr = 0,
  chrom_length = potato_chrom_length
)

## ----hidecan-nrows, fig.width = 10, fig.height = 8----------------------------
## Specifying the number of rows
hidecan_plot(
  gwas_list = x[["GWAS"]],          
  de_list = x[["DE"]],              
  can_list = x[["CAN"]],            
  score_thr_gwas = -log10(0.0001),
  score_thr_de = -log10(0.05),
  log2fc_thr = 0,
  n_rows = 3
)

## ----hidecan-ncols, fig.width = 10, fig.height = 10---------------------------
## Specifying the number of columns
hidecan_plot(
  gwas_list = x[["GWAS"]],          
  de_list = x[["DE"]],              
  can_list = x[["CAN"]],            
  score_thr_gwas = -log10(0.0001),
  score_thr_de = -log10(0.005),
  log2fc_thr = 0,
  n_cols = 3
)

## ----show-window-size-error, echo = FALSE, error = TRUE, fig.width = 1, fig.height = 1----
hidecan_plot(gwas_list = x[["GWAS"]],
             de_list = x[["DE"]],
             can_list = x[["CAN"]],
             score_thr_gwas = -log10(0.0001),
             score_thr_de = -log10(0.05),
             log2fc_thr = 0,
             label_size = 2)

## ----use-ggsave, eval = FALSE-------------------------------------------------
#  p <- hidecan_plot(
#    gwas_list = x[["GWAS"]],
#    de_list = x[["DE"]],
#    can_list = x[["CAN"]],
#    score_thr_gwas = -log10(0.0001),
#    score_thr_de = -log10(0.05),
#    log2fc_thr = 0,
#    label_size = 2
#  )
#  
#  ggplot2::ggsave("hidecan_plot.pdf", p, width = 10, height = 10)


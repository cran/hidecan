## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(tibble)
library(dplyr)
library(purrr)

## ----setup--------------------------------------------------------------------
library(hidecan)
library(tibble)
library(dplyr)
library(purrr)

## ----get-example-data---------------------------------------------------------
x <- get_example_data()

str(x)

## ----gwas-data-constructor----------------------------------------------------
gwas_data <- GWAS_data(x[["GWAS"]])

class(gwas_data)

head(gwas_data)

## ----de-data-constructor------------------------------------------------------
de_data <- DE_data(x[["DE"]])

class(de_data)

head(de_data)

## ----can-data-constructor-----------------------------------------------------
## CAN_data constructor
can_data <- CAN_data(x[["CAN"]])

class(can_data)

head(can_data)

## ----gwas-data-no-chromosome-column, error = TRUE-----------------------------
gwas_wrong_input <- x[["GWAS"]] |> 
  select(-chromosome)

GWAS_data(gwas_wrong_input)

## ----show-de-data-new-columns-------------------------------------------------
## Input tibble
head(x[["DE"]])

## Output of the DE_data constructor
head(de_data)

## ----combine-chrom-length-----------------------------------------------------
chrom_length <- combine_chrom_length(list(gwas_data,
                                          de_data,
                                          can_data))

chrom_length

## ----compute-chrom-length-----------------------------------------------------
head(compute_chrom_length(gwas_data), 3)

head(compute_chrom_length(de_data), 3)

## ----apply-threshold-gwas-----------------------------------------------------
dim(gwas_data)

gwas_data_thr <- apply_threshold(gwas_data, 
                                 score_thr = 4)

class(gwas_data_thr)

dim(gwas_data_thr)

head(gwas_data_thr)

## ----apply-threshold-de-------------------------------------------------------
dim(de_data)

de_data_thr <- apply_threshold(de_data, 
                               score_thr = 2,
                               log2fc_thr = 0.5)

class(de_data_thr)

dim(de_data_thr)

head(de_data_thr)

## ----apply-threshold-can------------------------------------------------------
dim(can_data)

can_data_thr <- apply_threshold(can_data, 
                                score_thr = 2,
                                log2fc_thr = 0.5)

class(can_data_thr)

dim(can_data_thr)

head(can_data_thr)

## ----create-hidecan-plot, fig.width = 10, fig.height = 10---------------------
create_hidecan_plot(
  list(gwas_data_thr,
       de_data_thr,
       can_data_thr),
  chrom_length
)


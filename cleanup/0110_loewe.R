#housekeeping
rm(list=ls(all=TRUE))

setwd('/home/bulat/NCI/cleanup')

#required libraries
library("dplyr")
library('drc')
library('nleqslv')
library('readr')
library('reshape2')
library('tibble')

source('ReshapeForDB.R')
file <- '2609_reshaped.Loewe'
reshaped <- readRDS(file)

to.bind <- ReshapeForDB(reshaped)
calculated.synergy.rbound.datalist = do.call(rbind, to.bind)

saveRDS(object = calculated.synergy.rbound.datalist, file = '0110_rbound.Loewe')
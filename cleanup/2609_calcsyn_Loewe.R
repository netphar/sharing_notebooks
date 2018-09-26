#housekeeping
rm(list=ls(all=TRUE))

#required libraries
library("dplyr")
library('drc')
library('nleqslv')
library('readr')
library('reshape2')
library('tibble')

#required functions
source('correctingConc.R')
source('rowInsert.R')
source('ReshapeData.R')
source('exchange.R')
source('ZIP.R')
source('HSA.R')
source('Bliss.R')
source('Loewe.R')
source('BaselineCorrectionSD2.R')
source('FittingSingleDrug.R')
source('CalculateSynergy.R')

setwd('/home/bulat/NCI/cleanup')
file <- '2609_reshaped'
reshaped <- readRDS(file)

CalculateSynergy(reshaped, method = 'Loewe', correction = T, Emin = 0, Emax = NA) -> reshaped.Loewe
saveRDS(object = reshaped.Loewe, file = '2609_reshaped.Loewe')
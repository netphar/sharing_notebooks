#working dir server
setwd('/home/bulat/NCI/cleanup')


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
source('ReshapeForDB.R')
source('ggplotRegression.R')

file <- '2609_reshaped'
reshaped <- readRDS(file)

CalculateSynergy(temp.reshaped, method = 'ZIP', correction = T, Emin = 0, Emax = NA) -> temp.reshaped.ZIP


# saving after reshape
saveRDS(object = temp.reshaped.ZIP, file = paste0(format(Sys.Date(), "%m_%d_%Y_"), "reshaped_ZIP"))

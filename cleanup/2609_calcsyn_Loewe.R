#working dir server
setwd('/home/bulat/NCI/cleanup')


#housekeeping
rm(list=ls(all=TRUE))
options(show.error.messages = F)

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

file <- '10_06_2018_reshaped'
temp.reshaped <- readRDS(file)

CalculateSynergy(temp.reshaped, method = 'Loewe', correction = T, Emin = 0, Emax = NA) -> temp.reshaped.Loewe


# saving after reshape
saveRDS(object = temp.reshaped.Loewe, file = paste0(format(Sys.Date(), "%m_%d_%Y_"), "reshaped_Loewe"))

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

file <- '10_06_2018_reshaped'
temp.reshaped <- readRDS(file)


options(show.error.messages = F)
CalculateSynergy(temp.reshaped, method = 'Bliss', correction = T, Emin = 0) -> temp.reshaped.Bliss

saveRDS(object = temp.reshaped.Bliss, file = paste0(format(Sys.Date(), "%m_%d_%Y_"), "reshaped_Bliss"))

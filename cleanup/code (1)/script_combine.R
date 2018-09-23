library(openxlsx)
library(plyr)
library(reshape2)
library(drc)
setwd("C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\database\\")
rm(list=ls(all=TRUE))

drug = read.xlsx("drug.xlsx") # 8 columns, no synonym and iupac
cell = read.xlsx("cell_line.xlsx") # 6 columns

# combine response data, 19 columns
response_oneil = read.csv("C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\Oneil\\response_oneil_final.csv", stringsAsFactors = F)
response_licci = read.csv("C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\Licci\\response_licci_final.csv", stringsAsFactors = F)
response_all = rbind(response_oneil, response_licci)
response_all$response = round(response_all$response, 2) # keep only two digits
response_all$synergy_bliss = round(response_all$synergy_bliss, 2)
response_all$synergy_zip = round(response_all$synergy_zip, 2)
response_all$synergy_loewe = round(response_all$synergy_loewe, 2)
response_all$synergy_hsa = round(response_all$synergy_hsa, 2)

response_all$drug_row_id = drug$drug_id[match(response_all$drug_row, drug$drug_name)]
response_all$drug_col_id = drug$drug_id[match(response_all$drug_col, drug$drug_name)]
response_all$cell_line_id = cell$cell_line_id[match(response_all$cell_line_name, cell$cell_line_name)]

# response_nci = read.csv("C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\NCI\\response_nci_final.csv")
# start = max(response_licci$block_id, response_oneil$block_id)
# response_nci$block_id = start + response_nci$block_id
# response_all = rbind(response_oneil, response_licci, response_nci)

# quality check
identical(max(response_all$block_id), length(unique(response_all$block_id)))

# combine summary data
summary_oneil = read.csv("C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\Oneil\\summary_oneil.csv", stringsAsFactors = F)
summary_licci = read.csv("C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\Licci\\summary_licci.csv", stringsAsFactors = F)
summary_all = rbind(summary_oneil, summary_licci)

summary_all$css = round(as.numeric(summary_all$css), 2)
summary_all$dss_row = round(as.numeric(summary_all$dss_row), 2)
summary_all$dss_col = round(as.numeric(summary_all$dss_col), 2)
summary_all$css_row = round(as.numeric(summary_all$css_row), 2)
summary_all$css_col = round(as.numeric(summary_all$css_col), 2)
summary_all$synergy_zip = round(as.numeric(summary_all$synergy_zip), 2)
summary_all$synergy_bliss = round(as.numeric(summary_all$synergy_bliss), 2)
summary_all$synergy_loewe = round(as.numeric(summary_all$synergy_loewe), 2)
summary_all$synergy_hsa = round(as.numeric(summary_all$synergy_hsa), 2)
summary_all$drug_row_id = drug$drug_id[match(summary_all$drug_row, drug$drug_name)]
summary_all$drug_col_id = drug$drug_id[match(summary_all$drug_col, drug$drug_name)]
summary_all$cell_line_id = cell$cell_line_id[match(summary_all$cell_line_name, cell$cell_line_name)]

identical(max(summary_all$block_id), length(unique(summary_all$block_id)))

# combine curve data
curve_oneil = read.csv("C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\Oneil\\curve_oneil.csv", stringsAsFactors = F)
curve_licci = read.csv("C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\Licci\\curve_licci.csv", stringsAsFactors = F)
curve_all = rbind(curve_oneil, curve_licci)
curve_all$drug_row_id = drug$drug_id[match(curve_all$drug_row, drug$drug_name)]
curve_all$drug_col_id = drug$drug_id[match(curve_all$drug_col, drug$drug_name)]

identical(max(curve_all$block_id), length(unique(curve_all$block_id)))

# combine surface data
surface_oneil = read.csv("C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\Oneil\\surface_oneil.csv", stringsAsFactors = F)
surface_licci = read.csv("C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\Licci\\surface_licci.csv", stringsAsFactors = F)
surface_all = rbind(surface_oneil, surface_licci)
surface_all$drug_row_id = drug$drug_id[match(surface_all$drug_row, drug$drug_name)]
surface_all$drug_col_id = drug$drug_id[match(surface_all$drug_col, drug$drug_name)]

# quality check
identical(max(surface_all$block_id), length(unique(surface_all$block_id)))

# quality check
length(setdiff(unique(response_all$block_id), unique(curve_all$block_id))) == 0
length(setdiff(unique(curve_all$block_id), unique(surface_all$block_id))) == 0

# write to database
head(response_all)
write.csv(response_all[, c("block_id","conc_r","conc_c","response","synergy_zip","synergy_bliss","synergy_loewe","synergy_hsa","cell_line_id","drug_row_id","drug_col_id")], "C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\database\\response.csv", row.names = F)

head(summary_all)
write.csv(summary_all[, -which(names(summary_all) %in% c("drug_row", "drug_col", "cell_line_name"))], "C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\database\\summary.csv", row.names = F)

head(curve_all)
write.csv(curve_all[, -which(names(curve_all) %in% c("drug_row", "drug_col", "conc_r_unit","conc_c_unit"))], "C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\database\\curve.csv", row.names = F)

head(surface_all)
write.csv(surface_all[, -which(names(surface_all) %in% c("drug_row", "drug_col", "conc_r_unit","conc_c_unit"))], "C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\database\\surface.csv", row.names = F)

# write the surface.csv into three parts
write.csv(surface_all[c(1:1500000), -which(names(surface_all) %in% c("drug_row", "drug_col", "conc_r_unit", "conc_c_unit", "drug_row_id", "drug_col_id"))], "C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\database\\surface_part1.csv", row.names = F)
write.csv(surface_all[c(1500001:3000000), -which(names(surface_all) %in% c("drug_row", "drug_col", "conc_r_unit","conc_c_unit","drug_row_id", "drug_col_id"))], "C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\database\\surface_part2.csv", row.names = F)
write.csv(surface_all[c(3000001:nrow(surface_all)), -which(names(surface_all) %in% c("drug_row", "drug_col", "conc_r_unit","conc_c_unit", "drug_row_id", "drug_col_id"))], "C:\\Users\\Localadmin_jtang\\Dropbox\\drugcomb\\data\\database\\surface_part3.csv", row.names = F)


save.image("C:/Users/Localadmin_jtang/Dropbox/drugcomb/data/database/result.RData")


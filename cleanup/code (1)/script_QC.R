library(gplots)
library(lattice)
source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/smoothing.R')
source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/own_rank.R')
source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/ReshapeData2.R')
source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/FittingSingleDrug2.R')
source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/CalculateSynergy2.R')
source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/ZIP2.R')
source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/Loewe2.R')
source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/Bliss2.R')
source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/HSA2.R')
source('C:/Users/Localadmin_jtang/Dropbox/drugcomb/code/BaselineCorrectionSD2.R')
options(scipen = 999) # suppress the scientific notation

load("C:/Users/Localadmin_jtang/Dropbox/drugcomb/data/database/result.RData")
response_all$drug_col = as.character(response_all$drug_col)
response_all$drug_row = as.character(response_all$drug_row)
response_all$conc_c_unit = as.character(response_all$conc_c_unit)
response_all$conc_r_unit = as.character(response_all$conc_r_unit)

i = which(summary_all$drug_row=='MK-4827' & summary_all$drug_col=='GELDANAMYCIN' & summary_all$cell_line_name=='SW837')
i = 1288

# tmp=surface_all[which(surface_all$block_id==i),]
data = ReshapeData2(response_all[which(response_all$block_id==i),], data.type = "inhibition")
data

# try some NA values
data$dose.response.mats[[1]][1,1] = NA
data$dose.response.mats[[1]][1,2] = NA
data$dose.response.mats[[1]][2,1] = NA

# NA values not allowed
bliss = CalculateSynergy2(data, method = "Bliss", correction = T, nan.handle = "L4", Emin = 0, Emax = NA)
zip = CalculateSynergy2(data, method = "ZIP", correction = T, nan.handle = "L4", Emin = 0, Emax = NA)
loewe = CalculateSynergy2(data, method = "LOEWE", correction = T, nan.handle = "L4", Emin = 0, Emax = NA)
hsa = CalculateSynergy2(data, method = "HSA", correction = T, nan.handle = "L4", Emin = 0, Emax = NA)

summary_all[which(summary_all$block_id==i),]
sum(bliss$scores[[1]])/28
sum(zip$scores[[1]])/28
sum(loewe$scores[[1]])/28
sum(hsa$scores[[1]])/28

mat1 = smoothing(data$dose.response.mats[[1]])
my_palette <- colorRampPalette(c("blue", "black", "red"))
graphics.off()
heatmap.2(mat1, Rowv = F, Colv = F, dendrogram = 'none',trace = 'none',col = my_palette, density.info = "none")

# contour plot
levels <- seq(0, 70, by = 5)
col1 <- colorRampPalette(c("green", "#FFFFFF"))(length(which(levels <= 
                                                               0)))
col2 <- colorRampPalette(c("#FFFFFF", "red"))(length(which(levels >=                                                           0)))
col <- c(col1, col2[-1])
x.2D <- (1:dim(mat1)[1] - 1)/(dim(mat1)[1] - 1)
y.2D <- (1:dim(mat1)[2] - 1)/(dim(mat1)[2] - 1)
plot.new()
plot.window(asp = NA, xlim = range(x.2D), ylim = range(y.2D), 
            "", xaxs = "i", yaxs = "i")
.filled.contour(x.2D, y.2D, z = t(mat1), levels, 
                col = col)
box()

pdf("168.pdf")
plot(wireframe(t(mat1)))
dev.off()


# make sure that surface is unique in the primary index
index = paste(surface_all$conc_c, surface_all$conc_r, surface_all$block_id, sep=" ")
index_duplidated = duplicated(index)
length(index_duplidated) == 0
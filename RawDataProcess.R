library(spectrolab)
library(tidyverse)
library(readxl)
library(DescTools)
library(cowplot)
library(gridGraphics)

# Pre processing of Raw data
source("~/Desktop/Nicotiana/Final/Code/Functions.R")

# Set global parameters for base R plot and theme for ggplot
par(cex.main = 1.5,  # title size
    font.main = 1,   # title font (2 = bold)
    cex.lab = 1.4,   # label size
    font.lab = 1,    # label font (2 = bold)
    cex.axis = 1.4,  # axis text size
    font.axis = 1,
    mar = c(4,4,4,3) + 0.1,
    bg = "transparent")   

my_theme <- theme_classic()+
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        text=element_text(size=17),
        legend.text = element_text(size = 12),
        legend.background = element_rect(fill="transparent"))


####Arizona#####
setwd("~/Desktop/Field_AZ/")

files = Sys.glob('2019-07-[0-9]*/*.asd')
spec_Arizona = read_spectra(path=files, type = "target_reflectance", format="asd")

# Metadata
metadata_table_az = read_excel("Metadata_Field_AZ.xlsx", sheet = "Samples") %>%
  drop_na(File_end_number_start)

summary_az = tibble(filename=unlist(strsplit(files,"/"))[seq(2,length(files)*2,2)],
                 type = rep(c("WR","WRL","BR","BRL"),each=5,times=nrow(metadata_table_az)),
                 Genotype_ID = rep(metadata_table_az$Genotype_ID, each=20),
                 Plant_ID = rep(metadata_table_az$Plant_ID, each=20))

# Exclude the "ExternalWR", only keep the leave measurements
df_az = cbind(summary_az, as.data.frame(spec_Arizona)) %>%
  filter(Genotype_ID!='ExternalWR')
metadata_table_az = filter(metadata_table_az,Genotype_ID!='ExternalWR')

# Dectect errors and outliers with plotting and LOF scores
# Plot WR scans
plot(as_spectra(filter(df_az,type == "WR")[,6:2156])) 
df = filter(df_az,type == "WR")[,6:2156]
filter(df_az,type == "WR")[apply(df, 1, mean) < 0.9, ]$Plant_ID # Plant 137, 815
# Make supplementary figure S2a
plot(as_spectra(filter(df_az,type == "WR")[,6:2156]),main="WR, Field_AZ",col = "#000000",xlab="Wavelength (nm)", ylab = "Reflectance",
     xaxs = "i",ylim=c(0,1.01))
plot(as_spectra(filter(df_az,type == "WR")[apply(df, 1, mean) < 0.9, 6:2156]), col="#CC79A7", add = TRUE)
abline(v=c(1000,1800), lty = 3)
plot_regions(as_spectra(filter(df_az,type == "WR")[,6:2156]),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add = TRUE)
legend(1800,0.9,legend=c("Scans","Outliers"),
       col=c("#000000","#CC79A7"), bg="transparent", 
       lty=1, lwd=3,cex=1,box.lty=0,seg.len = 2,x.intersp=0.5, y.intersp = 0.5)
Fig_S5a <- recordPlot()

df_az = filter(df_az, Plant_ID != 137, Plant_ID != 815) # remove the plants
metadata_table_az = filter(metadata_table_az,Plant_ID != 137, Plant_ID != 815)
# Plot BR scans
plot(as_spectra(filter(df_az,type == "BR")[,6:2156])) 
df = filter(df_az,type == "BR")[,6:2156]
filter(df_az,type == "BR")[apply(df, 1, mean) > 0.1, ]$Plant_ID  # Plant 165,209,813
df_az = filter(df_az, Plant_ID != 165, Plant_ID != 209, Plant_ID != 813) # remove the plants
metadata_table_az = filter(metadata_table_az,Plant_ID != 165, Plant_ID != 209, Plant_ID != 813)
# Plot WRL
plot(as_spectra(filter(df_az,type == "WRL")[,6:2156])) # OK
# Plot BRL
plot(as_spectra(filter(df_az,type == "BRL")[,6:2156])) # OK

# LOF outlier detector
df_az <- detect_outliers(df_az, 6,5) 

# Calculated Reflectance and Absolute uncertainty 
results <- calculate_CR_AU(df_az[,6:2156])
CR <- results$CR
AU <- results$AU

metadata_table_az = select(metadata_table_az,Genotype_ID, Plant_ID,Leaf_Num,batch,Day,Ctime)

A = cbind(metadata_table_az,AU) 
colnames(A) = c('Genotype_ID', "Plant_ID", "Leaf_Num", "Batch","Day","Ctime",350:2500)
A$Site = "Field_AZ"

R = cbind(metadata_table_az,CR) 
colnames(R) = c('Genotype_ID', "Plant_ID", "Leaf_Num", "Batch", "Day","Ctime", 350:2500)
R$Site = "Field_AZ"

# Take a look at all calculated reflectance
plot(as_spectra(R[,7:2157]),ylab = "Calculated Reflectance", main="Field_AZ",xlab = "Wavelength",xaxs = "i", yaxs = "i")

A = select(A,7:2157,1:6,2158)
R = select(R,7:2157,1:6,2158)

# Save processed dataset
saveRDS(A, "~/Desktop/Nicotiana/Final/Data_rds/A_AZ.rds")
saveRDS(R, "~/Desktop/Nicotiana/Final/Data_rds/R_AZ.rds")

####Jena#####
setwd("~/Desktop/Glasshouse")
files = Sys.glob('Schuman_20181213/*.asd')
spec_Jena = read_spectra(path=files, type = "target_reflectance", format="asd")

# Metadata
gene_id2 = read_excel("Metadata_Glasshouse.xlsx", sheet = "SampleID") %>%
  dplyr::select(Plant_ID, Genotype_ID) %>%
  filter(Genotype_ID != "irMAX2")

metadata_table_jena = read_excel("Metadata_Glasshouse.xlsx", sheet = "Samples") %>%
  dplyr::select(Plant_ID, batch,Ctime, File_end_number_start,WR_offset,'Leaf-WR_offset',BR_offset,'Leaf-BR_offset') %>%
  mutate(Plant_ID = as.numeric(Plant_ID)) %>%
  drop_na(Plant_ID) %>%
  left_join(.,gene_id2, by="Plant_ID")%>%
  drop_na(Plant_ID,Genotype_ID)

filename = unlist(strsplit(files,"/"))[seq(2,length(files)*2,2)]

summary_jena = tibble(filename=filename,
                 type = rep(c("WR","WRL","BR","BRL"),each=5,times=length(filename)/20),
                 Genotype_ID = rep(metadata_table_jena$Genotype_ID, each=20),
                 Plant_ID = rep(metadata_table_jena$Plant_ID, each=20),
                 batch = rep(metadata_table_jena$batch, each=20),
                 Ctime = rep(metadata_table_jena$Ctime, each=20))

df_jena = cbind(summary_jena, as.data.frame(spec_Jena))
metadata_jena = select(metadata_table_jena,Genotype_ID, Plant_ID,batch,Ctime)

# Dectect errors and outliers with plotting and LOF scores
# Plot WR scans
plot(as_spectra(filter(df_jena,type == "WR")[,8:2158])) # OK
# Plot BR scans
plot(as_spectra(filter(df_jena,type == "BR")[,8:2158])) # OK
# Plot WRL scans
plot(as_spectra(filter(df_jena,type == "WRL")[,8:2158])) # OK
# Plot BRL scans
plot(as_spectra(filter(df_jena,type == "BRL")[,8:2158])) # OK

# Make supplementary figure S2b
plot(as_spectra(filter(df_jena,type == "BR")[,8:2158]),main="BR, Glasshouse",col = "#000000",xlab="Wavelength (nm)", ylab = "Reflectance",
     xaxs = "i")
plot(as_spectra(filter(df_jena,filename %in% c("Jena_pre-tests_post-UVB_00073.asd","Jena_pre-tests_post-UVB_00554.asd"))[, 8:2158]), col="#CC79A7", add = TRUE)
abline(v=c(1000,1800), lty = 3)
plot_regions(as_spectra(df_jena[,8:2158]),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add = TRUE)
legend(1800,0.04,legend=c("Scans","Outliers"),
       col=c("#000000","#CC79A7"), bg="transparent", 
       lty=1, lwd=3,cex=1,box.lty=0,seg.len = 2,x.intersp=0.5, y.intersp = 0.5)
Fig_S5b <- recordPlot()

# LOF outlier detector
df_jena <- detect_outliers(df_jena, 8,5) #scan 403,073,554 set to NAs

# Calculated Reflectance and Absolute uncertainty 
results <- calculate_CR_AU(df_jena[,8:2158])
CR <- results$CR
AU <- results$AU

A_2 = cbind(metadata_jena,AU) 
colnames(A_2) = c('Genotype_ID', "Plant_ID", "Batch","Ctime",350:2500)
A_2$Site = "Glasshouse"

R_2= cbind(metadata_jena,CR) 
colnames(R_2) = c('Genotype_ID', "Plant_ID","Batch", "Ctime", 350:2500)
R_2$Site = "Glasshouse"

# Take a look at all calculated reflectance
plot(as_spectra(R_2[,5:2155]),ylab = "Calculated Reflectance", main = "Glasshouse", xlab = "Wavelength",xaxs = "i", yaxs = "i")

A_2 = select(A_2,5:2155,1:4,2156)
R_2 = select(R_2,5:2155,1:4,2156)

# Save processed dataset
saveRDS(A_2, "~/Desktop/Nicotiana/Final/Data_rds/A_Jena.rds")
saveRDS(R_2, "~/Desktop/Nicotiana/Final/Data_rds/R_Jena.rds")

####Utah#####
setwd("~/Desktop/Field_UT")
## Lineages
files_am = Sys.glob('2019-05-24/lineages/*.asd')
spec_Utah_am = read_spectra(path=files_am, type = "target_reflectance", format="asd") 

files_noon = Sys.glob('2019-05-27/lineages/*.asd')
spec_Utah_noon = read_spectra(path=files_noon, type = "target_reflectance", format="asd")

files_pm = Sys.glob('2019-05-25/lineages/*.asd')
spec_Utah_pm = read_spectra(path=files_pm, type = "target_reflectance", format="asd") 

gene_id3 = read_excel("Metadata_Field_UT.xlsx", sheet = "Lineages") %>%
  filter(Genotype_ID!="irMAX2.1",Genotype_ID!="irMAX2.2")

summary_ut = tibble(filename=unlist(strsplit(files_am,"/"))[seq(3,length(files_am)*3,3)],
                    type = rep(c("WR","WRL","BR","BRL"),each=5,times=nrow(gene_id3)),
                    Genotype_ID = rep(gene_id3$Genotype_ID, each=20),
                    Plant_ID = rep(gene_id3$Plant_ID, each=20))

# Initial check on extreme values
# am
ut_lin_am = cbind(summary_ut, as.data.frame(spec_Utah_am))
ut_lin_am$Ctime = "am"
# Dectect errors and outliers with plotting and LOF scores
# Plot WR scans
plot(as_spectra(filter(ut_lin_am,type == "WR")[,6:2156])) # OK
# Plot BR scans
plot(as_spectra(filter(ut_lin_am,type == "BR")[,6:2156])) # OK
# Plot WRL scans
plot(as_spectra(filter(ut_lin_am,type == "WRL")[,6:2156])) 
filter(ut_lin_am,type == "WRL")[filter(ut_lin_am,type == "WRL")[,156] > 0.5,]$filename #"ref_01065.asd"
# Plot BRL scans
plot(as_spectra(filter(ut_lin_am,type == "BRL")[,6:2156]))
filter(ut_lin_am,type == "BRL")[filter(ut_lin_am,type == "BRL")[,656] < 0.1,]$filename # "ref_01375.asd"

# LOF outlier detector
ut_lin_am <- detect_outliers(ut_lin_am, 6,5) 

# noon
ut_lin_noon = cbind(summary_ut, as.data.frame(spec_Utah_noon))
ut_lin_noon$Ctime = "noon"

plot(as_spectra(filter(ut_lin_noon,type == "WR")[,6:2156])) 
plot(as_spectra(filter(ut_lin_noon,type == "BR")[,6:2156])) 
plot(as_spectra(filter(ut_lin_noon,type == "WRL")[,6:2156])) 
filter(ut_lin_noon,type == "WRL")[filter(ut_lin_noon,type == "WRL")[,156] > 0.3,]$filename #1305,1365
plot(as_spectra(filter(ut_lin_noon,type == "BRL")[,6:2156])) 
filter(ut_lin_noon,type == "BRL")[filter(ut_lin_noon,type == "BRL")[,656] < 0.1,]$filename # "ref_01255.asd" 

# LOF outlier detector
ut_lin_noon <- detect_outliers(ut_lin_noon, 6,5) # Remove plant 47, scan 784,1361,1381

# pm
ut_lin_pm = cbind(summary_ut, as.data.frame(spec_Utah_pm))
ut_lin_pm$Ctime = "pm"
plot(as_spectra(filter(ut_lin_pm,type == "WR")[,6:2156]))
plot(as_spectra(filter(ut_lin_pm,type == "BR")[,6:2156]))
filter(ut_lin_pm,type == "BR")[filter(ut_lin_pm,type == "BR")[,1156] > 0.1,]$filename # "ref_01090.asd" "ref_01190.asd" "ref_01210.asd"
plot(as_spectra(filter(ut_lin_pm,type == "WRL")[,6:2156]))
filter(ut_lin_pm,type == "WRL")[filter(ut_lin_pm,type == "WRL")[,1656] > 0.2,]$Plant_ID # remove 61,69
ut_lin_pm = filter(ut_lin_pm, Plant_ID != 61,Plant_ID != 69)

plot(as_spectra(filter(ut_lin_pm,type == "BRL")[,6:2156]))
filter(ut_lin_pm,type == "BRL")[filter(ut_lin_pm,type == "BRL")[,656] < 0.1,]$filename # "ref_00895.asd" "ref_00935.asd"

# LOF outlier detector
ut_lin_pm <- detect_outliers(ut_lin_pm, 6,5) # removed plant 21,42,scan 984,976

# Combine ut_lin dataset together
df_ut_lin = rbind(ut_lin_am,ut_lin_noon,ut_lin_pm)

# Calculated Reflectance and Absolute uncertainty 
results <- calculate_CR_AU(df_ut_lin[,6:2156])
CR <- results$CR
AU <- results$AU

meta_ut_lin = dplyr::select(df_ut_lin, Genotype_ID,Plant_ID,Ctime)[seq(1,nrow(df_ut_lin),20),]

A_3 = cbind(meta_ut_lin,AU)
colnames(A_3) = c('Genotype_ID', "Plant_ID", "Ctime",350:2500)
A_3$Site = "Field_UT"

R_3 = cbind(meta_ut_lin,CR)
colnames(R_3) = c('Genotype_ID', "Plant_ID", "Ctime",350:2500)
R_3$Site = "Field_UT"

plot(as_spectra(R_3[,4:2154]), main="Utah",ylab = "Calculated Reflectance", xlab = "Wavelength",xaxs = "i", yaxs = "i")
R_3[R_3[,1154] > 0.15,1:3] # Outliers, need to be removed. AU need to be recalculated
which(R_3[,1154] > 0.15) # 13,17,77,132,144

# Make supplementary figure S2c
plot(as_spectra(R_3[,4:2154]),main="Calculated Reflectance, Field_UT",col = "#000000",xlab="Wavelength (nm)", ylab = "Reflectance",
     xaxs = "i")
plot(as_spectra(R_3[c(13,17,77,132,144), 4:2154]), col="#CC79A7", add = TRUE)
abline(v=c(1000,1800), lty = 3)
plot_regions(as_spectra(R_3[,4:2154]),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add = TRUE)
legend(1800,0.7,legend=c("Plants","Outliers"),
       col=c("#000000","#CC79A7"), bg="transparent", 
       lty=1, lwd=3,cex=1,box.lty=0,seg.len = 2,x.intersp=0.5, y.intersp = 0.5)
Fig_S5c <- recordPlot()

pdf("~/Desktop/Nicotiana/Final/Fig/SFig5.pdf",width = 18, height=6)
plot_grid(Fig_S5a,Fig_S5b,Fig_S5c, 
          nrow = 1,labels = c('(a)','(b)','(c)'))
dev.off()

ut_lin_am = filter(ut_lin_am, Plant_ID != 31,Plant_ID != 35)
ut_lin_noon = filter(ut_lin_noon, Plant_ID != 42)
ut_lin_pm = filter(ut_lin_pm, Plant_ID != 47,Plant_ID != 59)

# Combine ut_lin dataset together
df_ut_lin = rbind(ut_lin_am,ut_lin_noon,ut_lin_pm)

# Calculated Reflectance and Absolute uncertainty 
results <- calculate_CR_AU(df_ut_lin[,6:2156])
CR <- results$CR
AU <- results$AU

meta_ut_lin = dplyr::select(df_ut_lin, Genotype_ID,Plant_ID,Ctime)[seq(1,nrow(df_ut_lin),20),]

A_3 = cbind(meta_ut_lin,AU)
colnames(A_3) = c('Genotype_ID', "Plant_ID", "Ctime",350:2500)
A_3$Site = "Field_UT"

R_3 = cbind(meta_ut_lin,CR)
colnames(R_3) = c('Genotype_ID', "Plant_ID", "Ctime",350:2500)
R_3$Site = "Field_UT"

plot(as_spectra(R_3[,4:2154]), main="Utah",ylab = "Calculated Reflectance", xlab = "Wavelength",xaxs = "i", yaxs = "i")

A_3 = dplyr::select(A_3,4:2154,1:3,2155)
R_3 = dplyr::select(R_3,4:2154,1:3,2155)

saveRDS(A_3, "~/Desktop/Nicotiana/Final/Data_rds/A_UT_lin.rds")
saveRDS(R_3, "~/Desktop/Nicotiana/Final/Data_rds/R_UT_lin.rds")

# Development
files_am = Sys.glob('2019-05-25/development/*.asd')
spec_Utah_am = read_spectra(path=files_am, type = "target_reflectance", format="asd")

files_noon = Sys.glob('2019-05-26/development/*.asd')
spec_Utah_noon = read_spectra(path=files_noon, type = "target_reflectance", format="asd")

gene_id4 = read_excel("Metadata_Field_UT.xlsx", sheet = "Development") %>%
  filter(Genotype_ID!="irMAX2.1",Genotype_ID!="irMAX2.2")

summary_utah = tibble(filename=unlist(strsplit(files_am,"/"))[seq(3,length(files_am)*3,3)],
                      type = rep(c("WR","WRL","BR","BRL"),each=5,times=nrow(gene_id4)),
                      Genotype_ID = rep(gene_id4$Genotype_ID, each=20),
                      Plant_ID = rep(gene_id4$Plant_ID, each=20),
                      Leaf_ID = rep(gene_id4$Leaf_ID, each=20))

# Initial check on extreme values
# am
ut_dev_am = cbind(summary_utah, as.data.frame(spec_Utah_am))
ut_dev_am$Ctime = "am"
plot(as_spectra(filter(ut_dev_am,type == "WR")[,7:2157])) 
plot(as_spectra(filter(ut_dev_am,type == "BR")[,7:2157])) 
plot(as_spectra(filter(ut_dev_am,type == "WRL")[,7:2157])) 
plot(as_spectra(filter(ut_dev_am,type == "BRL")[,7:2157])) 
filter(ut_dev_am,type == "BRL")[filter(ut_dev_am,type == "BRL")[,608] < 0.1,]$filename #"ref_00795.asd"

df_dev_am = mutate(ut_dev_am,Plant_ID = paste0(Plant_ID,"_",Leaf_ID))
df_dev_am <- detect_outliers(df_dev_am, 7,3) #remove 13_1,scan 602,629,579,659,779,819
df_dev_am = mutate(df_dev_am,Plant_ID = sub("_.*$", "", Plant_ID))

# noon
ut_dev_noon = cbind(summary_utah, as.data.frame(spec_Utah_noon))
ut_dev_noon$Ctime = "noon"
plot(as_spectra(filter(ut_dev_noon,type == "WR")[,7:2157])) 
plot(as_spectra(filter(ut_dev_noon,type == "BR")[,7:2157])) 
plot(as_spectra(filter(ut_dev_noon,type == "WRL")[,7:2157])) 
plot(as_spectra(filter(ut_dev_noon,type == "BRL")[,7:2157])) 

df_dev_noon = mutate(ut_dev_noon,Plant_ID = paste0(Plant_ID,"_",Leaf_ID))
df_dev_noon <- detect_outliers(df_dev_noon, 7,3) #remove 11_2,13_1, scan 696
df_dev_noon = mutate(df_dev_noon,Plant_ID = sub("_.*$", "", Plant_ID))

df_ut_dev = rbind(df_dev_am,df_dev_noon)

# Calculated Reflectance and Absolute uncertainty 
results <- calculate_CR_AU(df_ut_dev[,7:2157])
CR <- results$CR
AU <- results$AU

meta_ut_dev = select(df_ut_dev, Genotype_ID,Plant_ID,Leaf_ID,Ctime)[seq(1,nrow(df_ut_dev),20),]

A_4 = cbind(meta_ut_dev,AU)
colnames(A_4) = c('Genotype_ID', "Plant_ID", "Leaf_ID","Ctime",350:2500)
A_4$Site = "Field_UT"

R_4 = cbind(meta_ut_dev,CR)
colnames(R_4) = c('Genotype_ID', "Plant_ID", "Leaf_ID","Ctime",350:2500)
R_4$Site = "Field_UT"

plot(as_spectra(R_4[,5:2155]), main="Field_UT (development)",ylab = "Calculated Reflectance", xlab = "Wavelength",xaxs = "i", yaxs = "i")

A_4 = select(A_4,5:2155,1:4,2156)
R_4 = select(R_4,5:2155,1:4,2156)

saveRDS(A_4, "~/Desktop/Nicotiana/Final/Data_rds/A_UT_dev.rds")
saveRDS(R_4, "~/Desktop/Nicotiana/Final/Data_rds/R_UT_dev.rds")


######### Isolation Forest
df = filter(df_az,type=="WR")[,6:2156]

iso_forest <- isolation.forest(df)
anomaly_scores <- predict(iso_forest, df, output_score = TRUE)
filenames = filter(df_az,type =="WR")[anomaly_scores > 0.8,]$filename # 52,65
filter(df_az,type =="WR")[anomaly_scores > 0.8,]$Plant_ID

###### z-score
z_scores <- scale(df)
# Outliers are typically defined as values with a z-score > 3 or < -3
outlier_zscore <- apply(z_scores, 1, function(x) any(abs(x) > 3))
# Filter out outliers
filter(ut_lin_pm,type == "BRL")[outlier_zscore, ]$filename # 

library(spectrolab)
library(tidyverse)
library(readxl)

# Pre processing of Raw data
####Arizona#####
setwd("~/Desktop/PhD/WCCER/Arizona")

files = Sys.glob('2019-07-[0-9]*/*.asd')
spec_Arizona = read_spectra(path=files, type = "target_reflectance", format="asd")

# Metadata
metadata_table_az = read_excel("Metadata_mature_plants_MCS.xlsx", sheet = "Samples") %>%
  drop_na(File_end_number_start)

summary_az = tibble(filename=unlist(strsplit(files,"/"))[seq(2,length(files)*2,2)],
                 type = rep(c("WR","WRL","BR","BRL"),each=5,times=nrow(metadata_table_az)),
                 Genotype_ID = rep(metadata_table_az$Genotype_ID, each=20),
                 Plant_ID = rep(metadata_table_az$Plant_ID, each=20))

# Exclude the "ExternalWR", only keep the leave measurements
df_az = cbind(summary_az, as.data.frame(spec_Arizona)) %>%
  filter(Genotype_ID!='ExternalWR')
metadata_table_az = filter(metadata_table_az,Genotype_ID!='ExternalWR')

# Initial check on extreme values
plot(as_spectra(filter(df_az,type == "WR")[,8:2158])) 
filter(df_az,type == "WR")[filter(df_az,type == "WR")[,708] < 0.9,]$Plant_ID # Plant 137, 815
df_az = filter(df_az, Plant_ID != 137, Plant_ID != 815) # remove the plants
metadata_table_az = filter(metadata_table_az,Plant_ID != 137, Plant_ID != 815)
plot(as_spectra(filter(df_az,type == "WR")[,8:2158])) 

plot(as_spectra(filter(df_az,type == "BR")[,8:2158])) 
filter(df_az,type == "BR")[filter(df_az,type == "BR")[,708] > 0.1,]$Plant_ID # Plant 165,209,813
df_az = filter(df_az, Plant_ID != 165, Plant_ID != 209, Plant_ID != 813) # remove the plants
metadata_table_az = filter(metadata_table_az,Plant_ID != 165, Plant_ID != 209, Plant_ID != 813)
plot(as_spectra(filter(df_az,type == "BR")[,8:2158])) 

plot(as_spectra(filter(df_az,type == "WRL")[,8:2158])) #looks fine
plot(as_spectra(filter(df_az,type == "BRL")[,8:2158])) #looks fine

# Calculated Reflectance and Absolute uncertainty 
AU = matrix(, nrow = nrow(df_az)/20, 2151)
CR = matrix(, nrow = nrow(df_az)/20, 2151)
n = 1
for (i in seq(1,nrow(df_az),20)){
  for (j in 6:ncol(df_az)){
    WR = mean(df_az[(i+1):(i+4),j])
    WRL = mean(df_az[(i+6):(i+9),j])
    BR = mean(df_az[(i+11):(i+14),j])
    BRL = mean(df_az[(i+16):(i+19),j])
    STD_WR = sd(df_az[(i+1):(i+4),j])
    STD_WRL = sd(df_az[(i+6):(i+9),j])
    STD_BR = sd(df_az[(i+11):(i+14),j])
    STD_BRL = sd(df_az[(i+16):(i+19),j])
    CR[n,j-5] = (BRL*WR - WRL*BR)/(WR-BR)
    AU[n,j-5] = sqrt((BR*(WRL-BRL)/(WR-BR)^2)^2 * (STD_WR/2)^2 + (BR/(WR-BR))^2 * (STD_WRL/2)^2 + (WR*(WRL-BRL)/(WR-BR)^2)^2 * (STD_BR/2)^2 + (WR/(WR-BR))^2 * (STD_BRL/2)^2)
  }
  n = n+1
}

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
saveRDS(A, "~/Desktop/Nicotiana/Final/A_AZ.rds")
saveRDS(R, "~/Desktop/Nicotiana/Final/R_AZ.rds")

####Jena#####
setwd("~/Desktop/PhD/WCCER/Jena")
files = Sys.glob('Schuman_20181213/*.asd')
spec_Jena = read_spectra(path=files, type = "target_reflectance", format="asd")
spec_Jena = spec_Jena[c(1:1299,1305:1325),]

# Metadata
gene_id2 = read_excel("Metadata_plants_20181213.xlsx", sheet = "SampleID") %>%
  dplyr::select(Plant_ID, Genotype_ID)

metadata_table_jena = read_excel("Metadata_plants_20181213.xlsx", sheet = "Samples") %>%
  dplyr::select(Plant_ID, batch,Ctime, File_end_number_start,WR_offset,'Leaf-WR_offset',BR_offset,'Leaf-BR_offset') %>%
  mutate(Plant_ID = as.numeric(Plant_ID)) %>%
  left_join(.,gene_id2, by="Plant_ID")%>%
  drop_na(File_end_number_start)

filename = unlist(strsplit(files,"/"))[seq(2,length(files)*2,2)]
filename = filename[c(1:1299,1305:1325)]

summary_jena = tibble(filename=filename,
                 type = rep(c("WR","WRL","BR","BRL"),each=5,times=length(filename)/20),
                 Genotype_ID = rep(metadata_table_jena$Genotype_ID, each=20),
                 Plant_ID = rep(metadata_table_jena$Plant_ID, each=20),
                 batch = rep(metadata_table_jena$batch, each=20),
                 Ctime = rep(metadata_table_jena$Ctime, each=20))

# Delete ones: two mismatch plants
df_jena = cbind(summary_jena, as.data.frame(spec_Jena))
df_jena = df_jena[!(is.na(summary_jena$Plant_ID)),]
metadata_jena = drop_na(metadata_table_jena,Plant_ID)%>%
  select(Genotype_ID, Plant_ID,batch,Ctime)

# Initial check on extreme values
plot(as_spectra(filter(df_jena,type == "WR")[,8:2158])) # looks fine
plot(as_spectra(filter(df_jena,type == "BR")[,8:2158])) # looks fine
plot(as_spectra(filter(df_jena,type == "WRL")[,8:2158])) #looks fine
plot(as_spectra(filter(df_jena,type == "BRL")[,8:2158])) 
filter(df_jena,type == "BRL")[filter(df_jena,type == "BRL")[,708] > 0.8,]$filename # ...01304.asd, last scan
filter(df_jena,type == "BRL")[filter(df_jena,type == "BRL")[,708] > 0.8,]$Plant_ID # Plant 90

# Calculated Reflectance and Absolute uncertainty 
AU = matrix(, nrow = nrow(df_jena)/20, 2151)
CR = matrix(, nrow = nrow(df_jena)/20, 2151)

n = 1
for (i in seq(1,nrow(df_jena),20)){
  for (j in 8:ncol(df_jena)){
    WR = mean(df_jena[(i+1):(i+4),j])
    WRL = mean(df_jena[(i+6):(i+9),j])
    BR = mean(df_jena[(i+11):(i+14),j])
    BRL = mean(df_jena[(i+16):(i+19),j])
    STD_WR = sd(df_jena[(i+1):(i+4),j])
    STD_WRL = sd(df_jena[(i+6):(i+9),j])
    STD_BR = sd(df_jena[(i+11):(i+14),j])
    STD_BRL = sd(df_jena[(i+16):(i+19),j])
    CR[n,j-7] = (BRL*WR - WRL*BR)/(WR-BR)
    AU[n,j-7] = sqrt((BR*(WRL-BRL)/(WR-BR)^2)^2 * (STD_WR/2)^2 + (BR/(WR-BR))^2 * (STD_WRL/2)^2 + (WR*(WRL-BRL)/(WR-BR)^2)^2 * (STD_BR/2)^2 + (WR/(WR-BR))^2 * (STD_BRL/2)^2)
  }
  n = n+1
}

# Recalculation for Plant 90 and remove the last scan of BRL
which(df_jena$Plant_ID==90) # 1241~1260
i = 1241
which(metadata_jena$Plant_ID==90) #63
n=63

for (j in 8:ncol(df_jena)){
  WR = mean(df_jena[(i+1):(i+4),j])
  WRL = mean(df_jena[(i+6):(i+9),j])
  BR = mean(df_jena[(i+11):(i+14),j])
  BRL = mean(df_jena[(i+16):(i+18),j]) # remove the last scan
  STD_WR = sd(df_jena[(i+1):(i+4),j])
  STD_WRL = sd(df_jena[(i+6):(i+9),j])
  STD_BR = sd(df_jena[(i+11):(i+14),j])
  STD_BRL = sd(df_jena[(i+16):(i+18),j]) # remove the last scan
  CR[n,j-7] = (BRL*WR - WRL*BR)/(WR-BR)
  AU[n,j-7] = sqrt((BR*(WRL-BRL)/(WR-BR)^2)^2 * (STD_WR/2)^2 + (BR/(WR-BR))^2 * (STD_WRL/2)^2 + (WR*(WRL-BRL)/(WR-BR)^2)^2 * (STD_BR/2)^2 + (WR/(WR-BR))^2 * (STD_BRL/2)^2)
}

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

# Remove irMAX2
A_2 = filter(A_2,Genotype_ID != "irMAX2")
R_2 = filter(R_2,Genotype_ID != "irMAX2")

# Save processed dataset
saveRDS(A_2, "~/Desktop/Nicotiana/Final/A_Jena.rds")
saveRDS(R_2, "~/Desktop/Nicotiana/Final/R_Jena.rds")

####Utah#####
setwd("~/Desktop/PhD/WCCER/Utah/UTAH_ASD18130_2019")
## Lineages
files_am = Sys.glob('2019-05-24/lineages/*.asd')
spec_Utah_am = read_spectra(path=files_am, type = "target_reflectance", format="asd") # Warning flag: ref_00060~64

files_noon = Sys.glob('2019-05-27/lineages/*.asd')
spec_Utah_noon = read_spectra(path=files_noon, type = "target_reflectance", format="asd")

files_pm = Sys.glob('2019-05-25/lineages/*.asd')
spec_Utah_pm = read_spectra(path=files_pm, type = "target_reflectance", format="asd") # Warning flag: ref_00840~41

gene_id3 = read_excel("../2019_Nicotiana_Utah_metadata_EC_MCS_CL.xlsx", sheet = "Lineages")

summary_ut = tibble(filename=unlist(strsplit(files_am,"/"))[seq(3,length(files_am)*3,3)],
                    type = rep(c("WR","WRL","BR","BRL"),each=5,times=nrow(gene_id3)),
                    Genotype_ID = rep(gene_id3$Genotype_ID, each=20),
                    Plant_ID = rep(gene_id3$Plant_ID, each=20))

# Initial check on extreme values
# am
ut_lin_am = cbind(summary_ut, as.data.frame(spec_Utah_am))
ut_lin_am$Ctime = "am"
plot(as_spectra(filter(ut_lin_am,type == "WR")[,6:2156])) 
filter(ut_lin_am,type == "WR")[filter(ut_lin_am,type == "WR")[,708] < 0.8,]$filename # ref_01000.asd first scan, will be removed so that's fine

plot(as_spectra(filter(ut_lin_am,type == "BR")[,6:2156])) 

plot(as_spectra(filter(ut_lin_am,type == "WRL")[,6:2156])) 
filter(ut_lin_am,type == "WRL")[filter(ut_lin_am,type == "WRL")[,1208] > 0.4,]$filename #ref_01065.asd, first scan, fine

plot(as_spectra(filter(ut_lin_am,type == "BRL")[,6:2156])) 
filter(ut_lin_am,type == "BRL")[filter(ut_lin_am,type == "BRL")[,608] < 0.3,]$filename #"ref_00178.asd"&"ref_00179.asd" last 2 scans, "ref_01375.asd" fine
filter(ut_lin_am,type == "BRL")[filter(ut_lin_am,type == "BRL")[,608] < 0.3,]$Plant_ID # remove 9
ut_lin_am = filter(ut_lin_am, Plant_ID != 9)

# noon
ut_lin_noon = cbind(summary_ut, as.data.frame(spec_Utah_noon))
ut_lin_noon$Ctime = "noon"
plot(as_spectra(filter(ut_lin_noon,type == "WR")[,6:2156])) 
filter(ut_lin_noon,type == "WR")[filter(ut_lin_noon,type == "WR")[,656] < 0.8,]$filename #"ref_00920.asd" "ref_00921.asd" "ref_00922.asd" "ref_00923.asd"
filter(ut_lin_noon,type == "WR")[filter(ut_lin_noon,type == "WR")[,656] < 0.8,]$Plant_ID # remove 47
ut_lin_noon = filter(ut_lin_noon, Plant_ID != 47)

plot(as_spectra(filter(ut_lin_noon,type == "BR")[,6:2156])) 
filter(ut_lin_noon,type == "BR")[filter(ut_lin_noon,type == "BR")[,656] > 0.1,]$filename #"ref_00470.asd" first fine "ref_00830.asd" first fine

plot(as_spectra(filter(ut_lin_noon,type == "WRL")[,6:2156])) 
filter(ut_lin_noon,type == "WRL")[filter(ut_lin_noon,type == "WRL")[,1156] > 0.2,]$filename #"ref_00065.asd" "ref_00066.asd" "ref_00067.asd" "ref_00068.asd" "ref_00069.asd" "ref_01305.asd" first fine"ref_01365.asd" first fine
filter(ut_lin_noon,type == "WRL")[filter(ut_lin_noon,type == "WRL")[,1156] > 0.2,]$Plant_ID #remove 4
filter(ut_lin_noon,type == "WRL")[filter(ut_lin_noon,type == "WRL")[,1906] > 0.15,]$filename #"ref_00205.asd" "ref_00286.asd" "ref_00287.asd" "ref_00288.asd" "ref_00289.asd" "ref_00825.asd" "ref_01305.asd" "ref_01365.asd"
filter(ut_lin_noon,type == "WRL")[filter(ut_lin_noon,type == "WRL")[,1906] > 0.15,]$Plant_ID # remove 15
ut_lin_noon = filter(ut_lin_noon, Plant_ID != 4, Plant_ID !=15)
plot(as_spectra(filter(ut_lin_noon,type == "WRL")[,6:2156])) 

plot(as_spectra(filter(ut_lin_noon,type == "BRL")[,6:2156])) 
filter(ut_lin_noon,type == "BRL")[filter(ut_lin_noon,type == "BRL")[,656] < 0.1,]$filename # "ref_01255.asd" fine
filter(ut_lin_noon,type == "BRL")[filter(ut_lin_noon,type == "BRL")[,1906] > 0.2,]$filename #"ref_00835.asd" "ref_00836.asd" "ref_00837.asd" "ref_00838.asd" "ref_00839.asd"
filter(ut_lin_noon,type == "BRL")[filter(ut_lin_noon,type == "BRL")[,1906] > 0.2,]$Plant_ID # remove 42
ut_lin_noon = filter(ut_lin_noon, Plant_ID != 42)
plot(as_spectra(filter(ut_lin_noon,type == "BRL")[,6:2156])) 

# pm
ut_lin_pm = cbind(summary_ut, as.data.frame(spec_Utah_pm))
ut_lin_pm$Ctime = "pm"
plot(as_spectra(filter(ut_lin_pm,type == "WR")[,6:2156]))
filter(ut_lin_pm,type == "WR")[filter(ut_lin_pm,type == "WR")[,1156] < 0.8,]$filename #"ref_00720.asd" "ref_01140.asd"

plot(as_spectra(filter(ut_lin_pm,type == "BR")[,6:2156]))
filter(ut_lin_pm,type == "BR")[filter(ut_lin_pm,type == "BR")[,1156] > 0.1,]$filename # "ref_01090.asd" "ref_01190.asd" "ref_01210.asd"
filter(ut_lin_pm,type == "BR")[which.min(filter(ut_lin_pm,type == "BR")[,1156]),]$filename #"ref_00830.asd"

plot(as_spectra(filter(ut_lin_pm,type == "WRL")[,6:2156]))
filter(ut_lin_pm,type == "WRL")[filter(ut_lin_pm,type == "WRL")[,1156] > 0.3,]$filename #"ref_00365.asd" "ref_00425.asd" "ref_00805.asd" "ref_00945.asd" "ref_01045.asd" "ref_01145.asd" "ref_01365.asd" "ref_01366.asd" "ref_01367.asd" "ref_01368.asd" "ref_01369.asd"
filter(ut_lin_pm,type == "WRL")[filter(ut_lin_pm,type == "WRL")[,1156] > 0.3,]$Plant_ID # remove 69
filter(ut_lin_pm,type == "WRL")[filter(ut_lin_pm,type == "WRL")[,1906] > 0.2,]$filename 
filter(ut_lin_pm,type == "WRL")[filter(ut_lin_pm,type == "WRL")[,1906] > 0.2,]$Plant_ID #remove 17, 61
ut_lin_pm = filter(ut_lin_pm, Plant_ID != 69, Plant_ID != 17, Plant_ID != 61)
plot(as_spectra(filter(ut_lin_pm,type == "WRL")[,6:2156]))

plot(as_spectra(filter(ut_lin_pm,type == "BRL")[,6:2156]))
filter(ut_lin_pm,type == "BRL")[filter(ut_lin_pm,type == "BRL")[,656] < 0.1,]$filename # "ref_00895.asd" "ref_00935.asd"
filter(ut_lin_pm,type == "BRL")[filter(ut_lin_pm,type == "BRL")[,1906] > 0.15,]$filename
filter(ut_lin_pm,type == "BRL")[filter(ut_lin_pm,type == "BRL")[,1906] > 0.15,]$Plant_ID # remove 1,2,47,59
ut_lin_pm = filter(ut_lin_pm, Plant_ID != 1, Plant_ID != 2, Plant_ID != 47, Plant_ID != 59)

df_ut_lin = rbind(ut_lin_am,ut_lin_noon,ut_lin_pm)

AU = matrix(, nrow = nrow(df_ut_lin)/20, 2151)
CR = matrix(, nrow = nrow(df_ut_lin)/20, 2151)

n = 1
for (i in seq(1,nrow(df_ut_lin),20)){
  for (j in 6:2156){
    WR = mean(df_ut_lin[(i+1):(i+4),j])
    WRL = mean(df_ut_lin[(i+6):(i+9),j])
    BR = mean(df_ut_lin[(i+11):(i+14),j])
    BRL = mean(df_ut_lin[(i+16):(i+19),j])
    STD_WR = sd(df_ut_lin[(i+1):(i+4),j])
    STD_WRL = sd(df_ut_lin[(i+6):(i+9),j])
    STD_BR = sd(df_ut_lin[(i+11):(i+14),j])
    STD_BRL = sd(df_ut_lin[(i+16):(i+19),j])
    CR[n,j-5] = (BRL*WR - WRL*BR)/(WR-BR)
    AU[n,j-5] = sqrt((BR*(WRL-BRL)/(WR-BR)^2)^2 * (STD_WR/2)^2 + (BR/(WR-BR))^2 * (STD_WRL/2)^2 + (WR*(WRL-BRL)/(WR-BR)^2)^2 * (STD_BR/2)^2 + (WR/(WR-BR))^2 * (STD_BRL/2)^2)
  }
  n = n+1
}

meta_ut_lin = dplyr::select(df_ut_lin, Genotype_ID,Plant_ID,Ctime)[seq(1,nrow(df_ut_lin),20),]

A_3 = cbind(meta_ut_lin,AU)
colnames(A_3) = c('Genotype_ID', "Plant_ID", "Ctime",350:2500)
A_3$Site = "Field_UT"

R_3 = cbind(meta_ut_lin,CR)
colnames(R_3) = c('Genotype_ID', "Plant_ID", "Ctime",350:2500)
R_3$Site = "Field_UT"

plot(as_spectra(R_3[,4:2154]), main="Utah",ylab = "Calculated Reflectance", xlab = "Wavelength",xaxs = "i", yaxs = "i")
which(R_3[,1156]>0.2) # 30,34
plot(as_spectra(R_3[30,4:2154]), col = "red", add=TRUE) #remove 31
plot(as_spectra(R_3[34,4:2154]), col = "red", add=TRUE) # remove 35
A_3 = A_3[-c(30,34),]
R_3 = R_3[-c(30,34),]

plot(as_spectra(R_3[,4:2154]), main="Field_UT (lineages)",ylab = "Calculated Reflectance", xlab = "Wavelength",xaxs = "i", yaxs = "i")

A_3 = dplyr::select(A_3,4:2154,1:3,2155)
R_3 = dplyr::select(R_3,4:2154,1:3,2155)

# Remove irMAX2
A_3 = filter(A_3,Genotype_ID!="irMAX2.1",Genotype_ID!="irMAX2.2")
R_3 = filter(R_3,Genotype_ID!="irMAX2.1",Genotype_ID!="irMAX2.2")

saveRDS(A_3, "~/Desktop/Nicotiana/Final/A_UT_lin.rds")
saveRDS(R_3, "~/Desktop/Nicotiana/Final/R_UT_lin.rds")

# Development
files_am = Sys.glob('2019-05-25/development/*.asd')
spec_Utah_am = read_spectra(path=files_am, type = "target_reflectance", format="asd")

files_noon = Sys.glob('2019-05-26/development/*.asd')
spec_Utah_noon = read_spectra(path=files_noon, type = "target_reflectance", format="asd")

gene_id4 = read_excel("../2019_Nicotiana_Utah_metadata_EC_MCS_CL.xlsx", sheet = "Development")

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
filter(ut_dev_am,type == "WR")[filter(ut_dev_am,type == "WR")[,708] < 0.85,]$filename # "ref_00340.asd" "ref_00560.asd" "ref_00700.asd"

plot(as_spectra(filter(ut_dev_am,type == "BR")[,7:2157])) 
filter(ut_dev_am,type == "BR")[filter(ut_dev_am,type == "BR")[,708] > 0.3,]$filename # "ref_00630.asd"

plot(as_spectra(filter(ut_dev_am,type == "WRL")[,7:2157])) 
filter(ut_dev_am,type == "WRL")[filter(ut_dev_am,type == "WRL")[,1208] > 0.4,]$filename #"ref_00185.asd" "ref_00245.asd"

plot(as_spectra(filter(ut_dev_am,type == "BRL")[,7:2157])) 
filter(ut_dev_am,type == "BRL")[filter(ut_dev_am,type == "BRL")[,608] < 0.1,]$filename #"ref_00795.asd"
filter(ut_dev_am,type == "BRL")[filter(ut_dev_am,type == "BRL")[,1907] > 0.15,]$filename #"ref_00135.asd" "ref_00136.asd" "ref_00137.asd" "ref_00138.asd" "ref_00139.asd" "ref_00159.asd" "ref_00735.asd" "ref_00736.asd" "ref_00737.asd"
filter(ut_dev_am,type == "BRL")[filter(ut_dev_am,type == "BRL")[,1907] > 0.15,]$Plant_ID # remove plant 3, 13
ut_dev_am = filter(ut_dev_am, Plant_ID != 3, Plant_ID !=13)

# noon
ut_dev_noon = cbind(summary_utah, as.data.frame(spec_Utah_noon))
ut_dev_noon$Ctime = "noon"
plot(as_spectra(filter(ut_dev_noon,type == "WR")[,7:2157])) 
filter(ut_dev_noon,type == "WR")[filter(ut_dev_noon,type == "WR")[,656] < 0.9,]$filename #"ref_00460.asd" "ref_00520.asd" "ref_00820.asd"

plot(as_spectra(filter(ut_dev_noon,type == "BR")[,7:2157])) 
filter(ut_dev_noon,type == "BR")[filter(ut_dev_noon,type == "BR")[,656] > 0.04,]$filename #"ref_00530.asd"

plot(as_spectra(filter(ut_dev_noon,type == "WRL")[,7:2157])) 
filter(ut_dev_noon,type == "WRL")[filter(ut_dev_noon,type == "WRL")[,1956] > 0.8,]$filename #"ref_00105.asd" "ref_00106.asd" "ref_00107.asd" "ref_00108.asd" "ref_00109.asd"
filter(ut_dev_noon,type == "WRL")[filter(ut_dev_noon,type == "WRL")[,1956] > 0.8,]$Plant_ID #remove 2
ut_dev_noon = filter(ut_dev_noon, Plant_ID != 2)
plot(as_spectra(filter(ut_dev_noon,type == "WRL")[,7:2157])) 

plot(as_spectra(filter(ut_dev_noon,type == "BRL")[,7:2157])) 

df_ut_dev = rbind(ut_dev_am,ut_dev_noon)

AU = matrix(, nrow = nrow(df_ut_dev)/20, 2151)
CR = matrix(, nrow = nrow(df_ut_dev)/20, 2151)
n = 1
for (i in seq(1,nrow(df_ut_dev),20)){
  for (j in 7:2157){
    WR = mean(df_ut_dev[(i+1):(i+4),j])
    WRL = mean(df_ut_dev[(i+6):(i+9),j])
    BR = mean(df_ut_dev[(i+11):(i+14),j])
    BRL = mean(df_ut_dev[(i+16):(i+19),j])
    STD_WR = sd(df_ut_dev[(i+1):(i+4),j])
    STD_WRL = sd(df_ut_dev[(i+6):(i+9),j])
    STD_BR = sd(df_ut_dev[(i+11):(i+14),j])
    STD_BRL = sd(df_ut_dev[(i+16):(i+19),j])
    CR[n,j-6] = (BRL*WR - WRL*BR)/(WR-BR)
    AU[n,j-6] = sqrt((BR*(WRL-BRL)/(WR-BR)^2)^2 * (STD_WR/2)^2 + (BR/(WR-BR))^2 * (STD_WRL/2)^2 + (WR*(WRL-BRL)/(WR-BR)^2)^2 * (STD_BR/2)^2 + (WR/(WR-BR))^2 * (STD_BRL/2)^2)
  }
  n = n+1
}

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

# Remove irMAX2
A_4 = readRDS("~/Desktop/Nicotiana/Final/A_UT_dev.rds")
R_4 = readRDS("~/Desktop/Nicotiana/Final/R_UT_dev.rds")
A_4 = filter(A_4,Genotype_ID!="irMAX2.1",Genotype_ID!="irMAX2.2")
R_4 = filter(R_4,Genotype_ID!="irMAX2.1",Genotype_ID!="irMAX2.2")

saveRDS(A_4, "~/Desktop/Nicotiana/Final/A_UT_dev.rds")
saveRDS(R_4, "~/Desktop/Nicotiana/Final/R_UT_dev.rds")

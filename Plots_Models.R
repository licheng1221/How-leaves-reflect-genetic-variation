library(MASS)
library(tidyverse)
library(spectrolab)
library(lme4)
library(nlme)
library(emmeans)
library(report)
library(cowplot)
library(gridGraphics)

# Load customized functions
source("~/Desktop/Nicotiana/Final/Code/Functions.R")
options(max.print = 3000)

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

setwd("~/Desktop/Nicotiana/Final/Fig")

### Read data and basic statistics
## Arizona
AZ = readRDS("~/Desktop/Nicotiana/Final/Data_rds/R_AZ.rds")
#AZ = readRDS("C:/Users/Meredith Schuman/Dropbox/Work with Merry/Nicotiana attenuata paper 1-Sep/Data/R_AZ.rds")

AZ_UT = filter(AZ, str_detect(Genotype_ID,"UT-WT"))
AZ_PL = filter(AZ,startsWith(Genotype_ID,"P"))
AZ_RIL = filter(AZ,startsWith(Genotype_ID,"M"))

## Jena
Jena =  readRDS("~/Desktop/Nicotiana/Final/Data_rds/R_Jena.rds")
#Jena =  readRDS("C:/Users/Meredith Schuman/Dropbox/Work with Merry/Nicotiana attenuata paper 1-Sep/Data/R_Jena.rds")

Jena[Jena$Genotype_ID=="AZ",]$Genotype_ID = "P_AZ"
Jena_Ref = filter(Jena, Genotype_ID %in% c("UT_WT","pRESC2NC","pSOL3NC"))
Jena_UT = filter(Jena, Genotype_ID %in% c("UT_WT"))
Jena_EV = filter(Jena, Genotype_ID %in% c("pRESC2NC","pSOL3NC"))
Jena_EV1= filter(Jena, Genotype_ID %in% c("pRESC2NC"))
Jena_EV2 = filter(Jena, Genotype_ID %in% c("pSOL3NC"))
Jena_PL = filter(Jena,startsWith(Genotype_ID,"P"))
Jena_TL = filter(Jena,startsWith(Genotype_ID,"ir"))

## Utah
UT_lin = readRDS("~/Desktop/Nicotiana/Final/Data_rds/R_UT_lin.rds")
UT_dev = readRDS("~/Desktop/Nicotiana/Final/Data_rds/R_UT_dev.rds")

#UT_lin = readRDS("C:/Users/Meredith Schuman/Dropbox/Work with Merry/Nicotiana attenuata paper 1-Sep/Data/R_UT_lin.rds")
#UT_dev = readRDS("C:/Users/Meredith Schuman/Dropbox/Work with Merry/Nicotiana attenuata paper 1-Sep/Data/R_UT_dev.rds")

# UT-WT from UT_dev, select Leaf 2 which is comparable with the other datasets
UT_UT = filter(UT_dev,Leaf_ID == 2) 
# EV from UT_lin
UT_EV = filter(UT_lin,Genotype_ID == "pRESC2NC")
# TLs from UT_lin
UT_TL = filter(UT_lin, Genotype_ID !="pRESC2NC")

## CV
CV_AZ_UT = cal_cv(AZ_UT[,51:2151])
CV_AZ_PL = cal_cv(AZ_PL[,51:2151])
CV_AZ_RIL = cal_cv(AZ_RIL[,51:2151])
CV_Jena_UT = cal_cv(Jena_UT[,51:2151])
CV_Jena_EV = cal_cv(Jena_EV[,51:2151])
CV_Jena_EV1 = cal_cv(Jena_EV1[,51:2151])
CV_Jena_EV2 = cal_cv(Jena_EV2[,51:2151])
CV_Jena_Ref = cal_cv(Jena_Ref[,51:2151])
CV_Jena_PL = cal_cv(Jena_PL[,51:2151])
CV_Jena_TL = cal_cv(Jena_TL[,51:2151])
CV_UT_EV = cal_cv(UT_EV[,51:2151])
CV_UT_UT = cal_cv(UT_UT[,51:2151])
CV_UT_TL = cal_cv(UT_TL[,51:2151])

### One-way ANOVA of UT-WT in three envs
All_UT = rbind(dplyr::select(AZ_UT,51:2151,"Site"),dplyr::select(UT_UT,51:2151,"Site"),dplyr::select(Jena_UT,51:2151,"Site"))

p = vector()
azut = vector()
azjena = vector()
utjena = vector()

for (i in 1:2101){
  res.aov <- aov(All_UT[,i]~All_UT$Site)
  p[i] = summary(res.aov)[[1]][["Pr(>F)"]][1]
  bc = pairwise.t.test(All_UT[,i],All_UT$Site, p.adjust.method = "bonferroni")
  pairwise = emmeans(res.aov, pairwise ~ Site)
  pairwise = as.data.frame(pairwise$contrasts %>% summary(infer=TRUE))
  azut[i] = pairwise$p.value[1]
  azjena[i] = pairwise$p.value[2]
  utjena[i] = pairwise$p.value[3]
}
# check the assumption: residual analysis
plot(res.aov) #OK

# Report model at 500 nm
report(aov(All_UT[,101]~All_UT$Site))

m1_p_o = t(as.data.frame(p))
colnames(m1_p_o) = 400:2500
m1_p_a = t(as.data.frame(p.adjust(p, method = "BY")))
colnames(m1_p_a) = 400:2500

azut_o = t(as.data.frame(azut))
colnames(azut_o) = 400:2500
azut_a =t(as.data.frame(p.adjust(azut, method = "BY")))
colnames(azut_a) = 400:2500

azjena_o = t(as.data.frame(azjena))
colnames(azjena_o) = 400:2500
azjena_a =t(as.data.frame(p.adjust(azjena, method = "BY")))
colnames(azjena_a) = 400:2500

utjena_o = t(as.data.frame(utjena))
colnames(utjena_o) = 400:2500
utjena_a =t(as.data.frame(p.adjust(utjena, method = "BY")))
colnames(utjena_a) = 400:2500

#plot(as_spectra(UT_UT[,51:2101]), col="#D55E00", main="ANOVA, UT-WT only",xlab="Wavelength (nm)", ylab = "Reflectance | CV | P-value",
#     ylim=c(0,1.01), xaxs = "i")
#plot(as_spectra(AZ_UT[,51:2101]), col="#56B4E9", add = TRUE)
#plot(as_spectra(Jena_UT[,51:2101]), col="#009E73", add = TRUE)
plot_quantile(as_spectra(UT_UT[,51:2101]), col= RGB("#D55E00"), main="ANOVA, UT-WT only",total_prob = 1, border=FALSE,xlab="Wavelength (nm)", ylab = "Reflectance | CV | P-value",
              ylim=c(0,1.01), xaxs = "i")
plot_quantile(as_spectra(AZ_UT[,51:2101]), col=RGB("#56B4E9"), total_prob = 1, border=FALSE,add = TRUE)
plot_quantile(as_spectra(Jena_UT[,51:2101]), col=RGB("#009E73"), total_prob = 1, border=FALSE,add = TRUE)
plot(as_spectra(CV_UT_UT), col="#D55E00",lwd=3, add = TRUE)
plot(as_spectra(CV_Jena_UT), col="#009E73", lwd=3,add = TRUE)
plot(as_spectra(CV_AZ_UT), col="#56B4E9", lwd=3,add = TRUE)
plot(as_spectra(m1_p_o),col = "#CC79A7",lty=3, lwd=2,add=TRUE)
plot(as_spectra(m1_p_a),col = "#CC79A7",lty=2,lwd=2,add=TRUE)
abline(v=c(1000,1800), lty = 3)
abline(h=0.05, lty = 3, col = "#CC79A7")
plot_regions(as_spectra(CV_UT_UT), col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add=TRUE)
legend(400,1,legend=c("Field_AZ","Field_UT","Glasshouse","CV F_AZ", "CV F_UT","CV Gh","P-value","Adj. p-value"),
       col=c(NA,NA,NA,"#56B4E9","#D55E00","#009E73","#CC79A7","#CC79A7"), bg="transparent", border="transparent",
       lty=c(NA,NA,NA,1,1,1,3,2), lwd = c(NA,NA,NA,3,3,3,2,2),fill=c(RGB("#56B4E9"),RGB("#D55E00"),RGB("#009E73"),NA,NA,NA,NA,NA),
       cex=1,box.lty=0,seg.len = 2,x.intersp=0.5, y.intersp=0.5)
Fig_4a <- recordPlot()

plot(as_spectra(azut_o), col="#000000", main="Post-hoc pairwise comparison",xlab="Wavelength (nm)", ylab = "P-value",
     xaxs = "i", lty=2, lwd=2,ylim=c(0,1.01))
plot(as_spectra(azjena_o), col="#56B4E9", lty=2,lwd=2,add = TRUE)
plot(as_spectra(utjena_o), col="#D55E00", lty=2,lwd=2,add = TRUE)
plot(as_spectra(azut_a), col="#000000", lwd=3,add = TRUE)
plot(as_spectra(azjena_a), col="#56B4E9", lwd=3,add = TRUE)
plot(as_spectra(utjena_a), col="#D55E00",lwd=3, add = TRUE)
abline(v=c(1000,1800), lty = 3)
abline(h=0.05, lty = 3, col = "#CC79A7")
plot_regions(as_spectra(p_o),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add = TRUE)
legend(1850,1,legend=c("AZ vs Gh","UT vs Gh","AZ vs UT"),
       col=c("#56B4E9","#D55E00","#000000"), bg="transparent", 
       lty=1, lwd=3, cex=1,box.lty=0,seg.len = 2,x.intersp=0.5,y.intersp = 0.5)
Fig_4b <- recordPlot()

pdf("Fig4.pdf",width = 12,height = 6)
plot_grid(Fig_4a,Fig_4b,nrow = 1,labels = c('(a)','(b)'))
dev.off()

### Comparison of Ref samples: small sample size, one-way ANOVA of EV1,EV2,UT-WT in Jena
p = vector()
for (i in 1:2101){
  res.aov <- aov(Jena_Ref[,i+50]~Genotype_ID, data=Jena_Ref)
  p[i] = summary(res.aov)[[1]][["Pr(>F)"]][1]
}

# check the assumption: residual analysis
plot(res.aov) #OK

# Report model at 500nm
report(aov(Jena_Ref[,151]~Genotype_ID, data=Jena_Ref))

m4_p_o = t(as.data.frame(p))
colnames(m4_p_o) = 400:2500
m4_p_a = t(as.data.frame(p.adjust(p, method = "BY")))
colnames(m4_p_a) = 400:2500

plot_quantile(as_spectra(Jena_UT[,51:2151]), col = RGB("#000000"), total_prob = 1, border=FALSE,main="UT-WT, EV1, EV2 in Glasshouse",ylab = "Reflectance | CV | P-value",
     ylim=c(0,1.01),xaxs = "i",xlab = "Wavelength (nm)")
plot_quantile(as_spectra(Jena_EV1[,51:2151]), col=RGB("#E69F00"),total_prob = 1, border=FALSE,add = TRUE)
plot_quantile(as_spectra(Jena_EV2[,51:2151]), col=RGB("#F0E442"), total_prob = 1, border=FALSE,add = TRUE)
plot(as_spectra(CV_Jena_UT), col="#000000", lwd=3,add = TRUE)
plot(as_spectra(CV_Jena_EV1), col="#E69F00",lwd=3, add = TRUE)
plot(as_spectra(CV_Jena_EV2), col="#F0E442",lwd=3, add = TRUE)
plot(as_spectra(m4_p_o), col="#CC79A7",lty=3,lwd=2,add = TRUE)
plot(as_spectra(m4_p_a), col="#CC79A7",lty=2,lwd=3,add = TRUE)
plot_regions(as_spectra(Jena_UT[,51:2151]), col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add=TRUE)
abline(v=c(1000,1800), lty = 3)
abline(h=0.05, lty = 3, col = "#CC79A7")
legend(1750,1.0,legend=c("UT-WT", "EV1","EV2","CV UT-WT","CV EV1","CV EV2","P-value","Adj. p-value"),
       col=c(NA,NA,NA,"#000000","#E69F00","#F0E442","#CC79A7","#CC79A7"), bg="transparent", border="transparent",
       lty=c(NA,NA,NA,1,1,1,3,2), lwd = c(NA,NA,NA,3,3,3,2,3),fill=c(RGB("#000000"),col=RGB("#E69F00"),col=RGB("#F0E442"),NA,NA,NA,NA,NA),
       cex=1,box.lty=0,seg.len = 2,x.intersp=0.5, y.intersp = 0.5)
Fig_S3a <- recordPlot()

### Linear regression
## AZ
az = mutate(AZ,GeneGroup = ifelse(Genotype_ID == "UT-WT", "aUT-WT",
                                  ifelse(startsWith(Genotype_ID,"P"),"PL","RIL"))) %>%
  mutate(Day = paste0("Day",Day),
         Leaf_Num = paste0("n_",Leaf_Num),
         Batch = paste0("batch_",Batch))

# Decide main effects on Wavelength 1877, 824, 693, 552
best_model_az(1877) # M1 < M2 < M3=M4
best_model_az(824)  # M3 = M4 < M2 < M1
best_model_az(693) # M3 = M4 < M2 < M1
best_model_az(552) # M3 = M4 < M2 < M1

# We select M3. Now check for interactions
step_interact_az(1877) 
step_interact_az(824) 
step_interact_az(693) 
step_interact_az(552) 
# No consistent results on adding interactions. So we stay with M3

# check the assumption: residual analysis
plot(lm(az[,1877-350+1] ~ GeneGroup + Ctime + Batch + Leaf_Num, data = az)) #OK

# Run linear regression across all wavelengths
p_genegroup = vector()
plril = vector()
plwt = vector()
rilwt = vector()
p_time = vector()
amnoon = vector()
ampm = vector()
noonpm = vector()
p_batch = vector()
p_leafnum = vector()

for (i in 1:2101){
  variable = az[,i+50]
  M = lm(variable ~ GeneGroup + Ctime + Batch + Leaf_Num, data = az) 
  A = anova(M)
  p_genegroup[i] = A$`Pr(>F)`[1]
  p_time[i] = A$`Pr(>F)`[2]
  p_batch[i] = A$`Pr(>F)`[3]
  p_leafnum[i] = A$`Pr(>F)`[4]
  pairwise = emmeans(M, pairwise ~ GeneGroup)
  pairwise = as.data.frame(pairwise$contrasts %>% summary(infer=TRUE))
  plwt[i] = pairwise$p.value[1]
  rilwt[i] = pairwise$p.value[2]
  plril[i] = pairwise$p.value[3]
  pairwise2 = emmeans(M, pairwise ~ Ctime)
  pairwise2 = as.data.frame(pairwise2$contrasts %>% summary(infer=TRUE))
  amnoon[i] = pairwise2$p.value[1]
  ampm[i] = pairwise2$p.value[2]
  noonpm[i] = pairwise2$p.value[3]
}

# Report model at 500nm
variable = az[,151]
report(lm(variable ~ GeneGroup + Ctime + Batch + Leaf_Num, data = az))

# Genotypic effects
m2_genegroup_o = t(as.data.frame(p_genegroup))
m2_genegroup_a = t(as.data.frame(p.adjust(p_genegroup, method = "BY")))
colnames(m2_genegroup_o) = 400:2500 
sigwavs(m2_genegroup_o) #400~1885
colnames(m2_genegroup_a) = 400:2500 
sigwavs(m2_genegroup_a) #409~760, 1001~1017,1128-1800

m2_plril_o = t(as.data.frame(plril))
colnames(m2_plril_o) = 400:2500
m2_plril_a = t(as.data.frame(p.adjust(plril, method = "BY")))
colnames(m2_plril_a) = 400:2500

m2_plwt_o = t(as.data.frame(plwt))
colnames(m2_plwt_o) = 400:2500
m2_plwt_a = t(as.data.frame(p.adjust(plwt, method = "BY")))
colnames(m2_plwt_a) = 400:2500
sigwavs(m2_plwt_a) #466-656

m2_rilwt_o = t(as.data.frame(rilwt))
colnames(m2_rilwt_o) = 400:2500
m2_rilwt_a = t(as.data.frame(p.adjust(rilwt, method = "BY")))
colnames(m2_rilwt_a) = 400:2500
sigwavs(m2_rilwt_a) #413-1378，1390-1528

plot_quantile(as_spectra(filter(az,GeneGroup == "RIL")[,51:2151]), col = RGB("#F0E442"), total_prob = 1, border=FALSE,main="Genotypic effects in Field_AZ",ylab = "Reflectance | P-value",
     ylim=c(0,1.01),xaxs = "i",xlab = "Wavelength (nm)")
plot_quantile(as_spectra(filter(az,GeneGroup == "aUT-WT")[,51:2151]), col=RGB("#000000"), total_prob = 1, border=FALSE,add = TRUE)
plot_quantile(as_spectra(filter(az,GeneGroup == "PL")[,51:2151]), col=RGB("#E69F00"),total_prob = 1, border=FALSE,add = TRUE)
plot(as_spectra(m2_genegroup_o), col="#CC79A7",lty=2,lwd=2,add = TRUE)
plot(as_spectra(m2_genegroup_a), col="#CC79A7",lwd=3,add = TRUE)
plot_regions(as_spectra(CV_AZ_RIL), col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add=TRUE)
abline(v=c(1000,1800), lty = 3)
abline(h=0.05, lty = 3, col = "#CC79A7")
legend(400,1.05,legend=c("PLs", "RILs","UT-WT","P-value","Adj. p-value"),
       col=c(NA,NA,NA,"#CC79A7","#CC79A7"), bg="transparent", border="transparent",
       lty=c(NA,NA,NA,2,1), lwd = c(NA,NA,NA,2,3),fill=c(RGB("#E69F00"),RGB("#F0E442"),col=RGB("#000000"),NA,NA),
       cex=1,box.lty=0,seg.len = 2,x.intersp=0.5, y.intersp = 0.5)
Fig_6a <-  recordPlot()

plot(as_spectra(m2_plril_o), col="#000000",main="Post-hoc pairwise comparision",xlab="Wavelength (nm)", ylab = "P-value",
     xaxs = "i",lty=2,lwd=2,ylim=c(0,1.01))
plot(as_spectra(m2_plwt_o), col="#E69F00", lty=2,lwd=2,add = TRUE)
plot(as_spectra(m2_rilwt_o), col="#F0E442",lty=2, lwd=2,add = TRUE)
plot(as_spectra(m2_plril_a), col="#000000", lwd=3,add = TRUE)
plot(as_spectra(m2_plwt_a), col="#E69F00", lwd=3,add = TRUE)
plot(as_spectra(m2_rilwt_a), col="#F0E442",lwd=3, add = TRUE)
abline(v=c(1000,1800), lty = 3)
abline(h=0.05, lty = 3, col = "#CC79A7")
plot_regions(as_spectra(CV_AZ_UT),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add = TRUE)
legend(1750,0.2,legend=c("PLs vs Ref","RILs vs Ref","PLs vs RILs"),
       col=c("#E69F00","#F0E442","#000000"), bg="transparent", 
       lty=1, lwd=3,cex=1,box.lty=0,seg.len = 2,x.intersp=0.5, y.intersp = 0.5)
Fig_6b <-  recordPlot()

# Comparison: UT-WT vs RIL in Arizona 
#Bootstrapping
n = 100
perm = sapply(1:n, function(i){
  boot_id = sample(nrow(AZ_RIL),nrow(AZ_UT),replace=FALSE)
  boot_AZ_RIL = AZ_RIL[boot_id,51:2151] 
  t(as.data.frame(apply(boot_AZ_RIL,2,sd) / apply(boot_AZ_RIL,2,mean)))
})

CV_boot_AZ_RIL = t(perm)
colnames(CV_boot_AZ_RIL) = 400:2500

plot_quantile(as_spectra(CV_boot_AZ_RIL), total_prob = 1, border = FALSE, 
              main="Bootstrapped CV, RILs vs Ref, Field_AZ",ylab = "CV",xlab="Wavelength (nm)",
              col = RGB("#F0E442"),xaxs = "i")
plot(mean(as_spectra(CV_boot_AZ_RIL)), col="#F0E442",lty=2,lwd=2,add = TRUE)
plot(as_spectra(CV_AZ_UT), col="#000000",lwd=3, add = TRUE)
abline(v=c(1000,1800), lty = 3)
plot_regions(as_spectra(CV_AZ_UT),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add = TRUE)
legend(1750,0.23,legend=c("RILs", "Mean RILs","Ref"),
       col=c(NA,"#F0E442","#000000"), bg="transparent", border="transparent",
       lty=c(NA,2,1), lwd = c(NA,2,3),fill=c(RGB("#F0E442"),NA,NA),
       cex=1,box.lty=0,seg.len = 2,x.intersp=0.5, y.intersp = 0.5)
Fig_6c <-  recordPlot()

# Comparison: UT-WT vs PL in Arizona 
#Bootstrapping
n = 100
perm = sapply(1:n, function(i){
  boot_id = sample(nrow(AZ_UT),nrow(AZ_PL),replace=FALSE)
  boot_AZ_UT = AZ_UT[boot_id,51:2151] 
  t(as.data.frame(apply(boot_AZ_UT,2,sd) / apply(boot_AZ_UT,2,mean)))
})

CV_boot_AZ_UT = t(perm)
colnames(CV_boot_AZ_UT) = 400:2500

plot_quantile(as_spectra(CV_boot_AZ_UT), total_prob = 1, border = FALSE, 
              main="Bootstrapped CV, PLs vs Ref, Field_AZ",ylab = "CV",xlab="Wavelength (nm)",
              col = RGB("#000000"), xaxs = "i")
plot(mean(as_spectra(CV_boot_AZ_UT)), col = "#000000",lty = 2,lwd=2,add = TRUE)
plot(as_spectra(CV_AZ_PL), col="#E69F00",lwd=3, add = TRUE)
abline(v=c(1000,1800), lty = 3)
plot_regions(as_spectra(CV_AZ_PL),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add = TRUE)
legend(1780,0.07,legend=c("Ref","mean Ref","PLs"),
       col=c(NA,"#000000","#E69F00"), bg="transparent", border="transparent",
       lty=c(NA,2,1), lwd = c(NA,2,3),fill=c(RGB("#000000"),NA,NA),
       cex=1,box.lty=0,seg.len = 2,x.intersp=0.5, y.intersp = 0.5)
Fig_6d <-  recordPlot()

# Time effect
m2_p_time_o = t(as.data.frame(p_time))
m2_p_time_a = t(as.data.frame(p.adjust(p_time, method = "BY")))
colnames(m2_p_time_o) = 400:2500 
colnames(m2_p_time_a) = 400:2500 
sigwavs(m2_p_time_a) #734-1350,1358-1420,1541-1694,1861-1904

plot_quantile(as_spectra(filter(az,Ctime=="am")[,51:2151]),ylim = c(0,1.01), total_prob = 1, border=FALSE,col = RGB("#F0E442"), main = "Effects of time in Field_AZ (M2)", 
              ylab = "Reflectance | P-value", xlab="Wavelength (nm)",xaxs = "i")
plot_quantile(as_spectra(filter(az,Ctime=="noon")[,51:2151]),col = RGB("#E69F00"), total_prob = 1, border=FALSE,add=TRUE)
plot_quantile(as_spectra(filter(az,Ctime=="pm")[,51:2151]),col = RGB("#0072B2"), total_prob = 1, border=FALSE,add=TRUE)
plot(as_spectra(m2_p_time_o),col = "#CC79A7", lwd=2, lty=2, add=TRUE)
plot(as_spectra(m2_p_time_a), col = "#CC79A7", lwd=3, add=TRUE)
abline(v=c(1000,1800), lty = 3)
abline(h=0.05, lty = 3, col="#CC79A7")
plot_regions(as_spectra(m2_p_time_o),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add=TRUE)
legend(1750,1,legend=c("Am", "Noon","Pm","P-value","Adj. p-value"),
       fill=c(RGB("#F0E442"),RGB("#E69F00"),RGB("#0072B2"),NA,NA), col=c(NA,NA,NA,"#CC79A7","#CC79A7"),bg="transparent", border="transparent",
       lty=c(NA,NA,NA,2,1), lwd=c(NA,NA,NA,2,3),cex=1,box.lty=0,seg.len = 2,x.intersp=0.5, y.intersp = 0.5)
Fig_8c <-  recordPlot()

m2_amnoon_o = t(as.data.frame(amnoon))
colnames(m2_amnoon_o) = 400:2500
m2_amnoon_a = t(as.data.frame(p.adjust(amnoon, method = "BY")))
colnames(m2_amnoon_a) = 400:2500

m2_ampm_o = t(as.data.frame(ampm))
colnames(m2_ampm_o) = 400:2500
m2_ampm_a = t(as.data.frame(p.adjust(ampm, method = "BY")))
colnames(m2_ampm_a) = 400:2500

m2_noonpm_o = t(as.data.frame(noonpm))
colnames(m2_noonpm_o) = 400:2500
m2_noonpm_a = t(as.data.frame(p.adjust(noonpm, method = "BY")))
colnames(m2_noonpm_a) = 400:2500

plot(as_spectra(m2_amnoon_o), col="#F0E442",main="Post-hoc, Effect of time in Field_AZ (M2)",xlab="Wavelength (nm)", ylab = "P-value",
     xaxs = "i",lty=2,lwd=2,ylim=c(0,1.01))
plot(as_spectra(m2_ampm_o), col="#0072B2", lty=2,lwd=2,add = TRUE)
plot(as_spectra(m2_noonpm_o), col="#E69F00",lty=2, lwd=2,add = TRUE)
plot(as_spectra(m2_amnoon_a), col="#F0E442", lwd=3,add = TRUE)
plot(as_spectra(m2_ampm_a), col="#0072B2", lwd=3,add = TRUE)
plot(as_spectra(m2_noonpm_a), col="#E69F00",lwd=3, add = TRUE)
abline(v=c(1000,1800), lty = 3)
abline(h=0.05, lty = 3, col = "#CC79A7")
plot_regions(as_spectra(CV_AZ_UT),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add = TRUE)
legend(400,0.3,legend=c("Am vs Noon","AM vs Pm","Noon vs Pm"),
       col=c("#F0E442","#0072B2","#E69F00"), bg="transparent", 
       lty=1, lwd=3,cex=1,box.lty=0,seg.len = 2,x.intersp=0.5, y.intersp = 0.5)
Fig_S4a <- recordPlot()

# Effect of leaf numbers
m2_leafnum_o = t(as.data.frame(p_leafnum))
m2_leafnum_a = t(as.data.frame(p.adjust(p_leafnum, method = "BY")))
colnames(m2_leafnum_o) = 400:2500
colnames(m2_leafnum_a) = 400:2500 
sigwavs(m2_leafnum_a) #508-641,693-707,731-1149,1400-1539,2344-2500

plot_quantile(as_spectra(filter(az,Leaf_Num=="n_1")[,51:2151]),ylim = c(0,1.01), total_prob = 1, border=FALSE,col = RGB("#E69F00"),main = "Effects of 1 vs 2 leaves in Field_AZ (M2)", 
              ylab = "Reflectance | P-value", xlab="Wavelength (nm)",xaxs = "i")
plot_quantile(as_spectra(filter(az,Leaf_Num=="n_2")[,51:2151]),col = RGB("#0072B2"), total_prob = 1, border=FALSE,add=TRUE)
plot(as_spectra(m2_leafnum_o),col = "#CC79A7", lwd=2, lty=2, add=TRUE)
plot(as_spectra(m2_leafnum_a), col = "#CC79A7", lwd=3, add=TRUE)
abline(v=c(1000,1800), lty = 3)
abline(h=0.05, lty = 3, col="#CC79A7")
plot_regions(as_spectra(m2_leafnum_o),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add=TRUE)
legend(700,1,legend=c("1 leaf","2 leaves","P-value","Adj. p-value"),
       col=c(NA,NA,"#CC79A7","#CC79A7"),fill=c(RGB("#E69F00"),RGB("#0072B2"),NA,NA), bg="transparent",border="transparent", 
       lty=c(NA,NA,2,1), lwd=c(NA,NA,2,3),cex=1,box.lty=0,seg.len = 2,x.intersp=0.5, y.intersp = 0.5)
Fig_9a <- recordPlot()

# Effect of batches
m2_batch_o = t(as.data.frame(p_batch))
m2_batch_a = t(as.data.frame(p.adjust(p_batch, method = "BY")))
colnames(m2_batch_o) = 400:2500
colnames(m2_batch_a) = 400:2500
sigwavs(m2_batch_a) #400-1897,2002-2500

plot(as_spectra(m2_batch_o), main = "Effects of batch in Field_AZ (M2)", ylim = c(0,0.1),col = "#CC79A7", lty=2, lwd=2, ylab = "P-value", xlab="Wavelength (nm)",xaxs = "i")
plot(as_spectra(m2_batch_a), lwd=3, col = "#CC79A7", add=TRUE)
abline(v=c(1000,1800), lty = 3)
abline(h=0.05, lty = 3, col="#CC79A7")
plot_regions(as_spectra(m2_batch_o),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add=TRUE)
legend(400,0.1,legend=c("P-value","Adj. p-value"),
       col=c("#CC79A7","#CC79A7"), bg="transparent", 
       lty=c(2,1), lwd=c(2,3),cex=1,box.lty=0,seg.len = 2,x.intersp=0.5, y.intersp = 0.5)
Fig_S4c <- recordPlot()

### Linear regression
## AZ_UT-WT
az_ut <- filter(az,Genotype_ID == "UT-WT")
# check the assumption: residual analysis
plot(lm(az_ut[,1877-350+1] ~ Ctime + Batch + Leaf_Num, data = az_ut)) #OK

# Run linear regression across all wavelengths
p_time = vector()
amnoon = vector()
ampm = vector()
noonpm = vector()
p_batch = vector()
p_leafnum = vector()

for (i in 1:2101){
  variable = az_ut[,i+50]
  M = lm(variable ~ Ctime + Batch + Leaf_Num, data = az_ut) 
  A = anova(M)
  p_time[i] = A$`Pr(>F)`[1]
  p_batch[i] = A$`Pr(>F)`[2]
  p_leafnum[i] = A$`Pr(>F)`[3]
  pairwise = emmeans(M, pairwise ~ Ctime)
  pairwise = as.data.frame(pairwise$contrasts %>% summary(infer=TRUE))
  amnoon[i] = pairwise$p.value[1]
  ampm[i] = pairwise$p.value[2]
  noonpm[i] = pairwise$p.value[3]
}

# Report model at 500nm
variable = az_ut[,151]
report(lm(variable ~ Ctime + Batch + Leaf_Num, data = az_ut))

# Time effect
m3_p_time_o = t(as.data.frame(p_time))
m3_p_time_a = t(as.data.frame(p.adjust(p_time, method = "BY")))
colnames(m3_p_time_o) = 400:2500 
colnames(m3_p_time_a) = 400:2500 
sigwavs(m3_p_time_a) #None

plot_quantile(as_spectra(filter(az_ut,Ctime=="am")[,51:2151]),ylim = c(0,1.01), total_prob = 1, border=FALSE,col = RGB("#F0E442"), main = "Effects of time in Field_AZ (M3)", 
              ylab = "Reflectance | P-value", xlab="Wavelength (nm)",xaxs = "i")
plot_quantile(as_spectra(filter(az_ut,Ctime=="noon")[,51:2151]),col = RGB("#E69F00"), total_prob = 1, border=FALSE,add=TRUE)
plot_quantile(as_spectra(filter(az_ut,Ctime=="pm")[,51:2151]),col = RGB("#0072B2"), total_prob = 1, border=FALSE,add=TRUE)
plot(as_spectra(m3_p_time_o),col = "#CC79A7", lwd=2, lty=2, add=TRUE)
plot(as_spectra(m3_p_time_a), col = "#CC79A7", lwd=3, add=TRUE)
abline(v=c(1000,1800), lty = 3)
abline(h=0.05, lty = 3, col="#CC79A7")
plot_regions(as_spectra(m3_p_time_o),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add=TRUE)
legend(1750,0.5,legend=c("Am", "Noon","Pm","P-value","Adj. p-value"),
       fill=c(RGB("#F0E442"),RGB("#E69F00"),RGB("#0072B2"),NA,NA), col=c(NA,NA,NA,"#CC79A7","#CC79A7"),bg="transparent", border="transparent",
       lty=c(NA,NA,NA,2,1), lwd=c(NA,NA,NA,2,3),cex=1,box.lty=0,seg.len = 2,x.intersp=0.5, y.intersp = 0.5)
Fig_8d <-  recordPlot()

m3_amnoon_o = t(as.data.frame(amnoon))
colnames(m3_amnoon_o) = 400:2500
m3_amnoon_a = t(as.data.frame(p.adjust(amnoon, method = "BY")))
colnames(m3_amnoon_a) = 400:2500

m3_ampm_o = t(as.data.frame(ampm))
colnames(m3_ampm_o) = 400:2500
m3_ampm_a = t(as.data.frame(p.adjust(ampm, method = "BY")))
colnames(m3_ampm_a) = 400:2500

m3_noonpm_o = t(as.data.frame(noonpm))
colnames(m3_noonpm_o) = 400:2500
m3_noonpm_a = t(as.data.frame(p.adjust(noonpm, method = "BY")))
colnames(m3_noonpm_a) = 400:2500

plot(as_spectra(m3_amnoon_o), col="#F0E442",main="Post-hoc, Effect of time in Field_AZ (M3)",xlab="Wavelength (nm)", ylab = "P-value",
     xaxs = "i",lty=2,lwd=2,ylim=c(0,1.01))
plot(as_spectra(m3_ampm_o), col="#0072B2", lty=2,lwd=2,add = TRUE)
plot(as_spectra(m3_noonpm_o), col="#E69F00",lty=2, lwd=2,add = TRUE)
plot(as_spectra(m3_amnoon_a), col="#F0E442", lwd=3,add = TRUE)
plot(as_spectra(m3_ampm_a), col="#0072B2", lwd=3,add = TRUE)
plot(as_spectra(m3_noonpm_a), col="#E69F00",lwd=3, add = TRUE)
abline(v=c(1000,1800), lty = 3)
abline(h=0.05, lty = 3, col = "#CC79A7")
plot_regions(as_spectra(m3_ampm_o),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add = TRUE)
legend(1750,0.3,legend=c("Am vs Noon","AM vs Pm","Noon vs Pm"),
       col=c("#F0E442","#0072B2","#E69F00"), bg="transparent", 
       lty=1, lwd=3,cex=1,box.lty=0,seg.len = 2,x.intersp=0.5, y.intersp = 0.5)
Fig_S4b <- recordPlot()

# Effect of leaf numbers
m3_leafnum_o = t(as.data.frame(p_leafnum))
m3_leafnum_a = t(as.data.frame(p.adjust(p_leafnum, method = "BY")))
colnames(m3_leafnum_o) = 400:2500
colnames(m3_leafnum_a) = 400:2500 

plot_quantile(as_spectra(filter(az_ut,Leaf_Num=="n_1")[,51:2151]),ylim = c(0,1.01), total_prob = 1, border=FALSE,col = RGB("#E69F00"),main = "Effects of 1 vs. 2 leaves in Field_AZ (M3)", 
              ylab = "Reflectance | P-value", xlab="Wavelength (nm)",xaxs = "i")
plot_quantile(as_spectra(filter(az_ut,Leaf_Num=="n_2")[,51:2151]),col = RGB("#0072B2"), total_prob = 1, border=FALSE,add=TRUE)
plot(as_spectra(m3_leafnum_o),col = "#CC79A7", lwd=2, lty=2, add=TRUE)
plot(as_spectra(m3_leafnum_a), col = "#CC79A7", lwd=3, add=TRUE)
abline(v=c(1000,1800), lty = 3)
abline(h=0.05, lty = 3, col="#CC79A7")
plot_regions(as_spectra(m3_leafnum_o),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add=TRUE)
legend(1750,0.45,legend=c("1 leaf","2 leaves","P-value","Adj. p-value"),
       col=c(NA,NA,"#CC79A7","#CC79A7"),fill=c(RGB("#E69F00"),RGB("#0072B2"),NA,NA), bg="transparent",border="transparent", 
       lty=c(NA,NA,2,1), lwd=c(NA,NA,2,3),cex=1,box.lty=0,seg.len = 2,x.intersp=0.5, y.intersp = 0.5)
Fig_9b <- recordPlot()

# Effect of batches
m3_batch_o = t(as.data.frame(p_batch))
m3_batch_a = t(as.data.frame(p.adjust(p_batch, method = "BY")))
colnames(m3_batch_o) = 400:2500
colnames(m3_batch_a) = 400:2500

plot(as_spectra(m3_batch_o), main = "Effects of batch in Field_AZ (M3)", ylim = c(0,1.01),col = "#CC79A7", lty=2, lwd=2, ylab = "P-value", xlab="Wavelength (nm)",xaxs = "i")
plot(as_spectra(m3_batch_a), lwd=3, col = "#CC79A7", add=TRUE)
abline(v=c(1000,1800), lty = 3)
abline(h=0.05, lty = 3, col="#CC79A7")
plot_regions(as_spectra(m3_batch_o),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add=TRUE)
legend(1800,1,legend=c("P-value","Adj. p-value"),
       col=c("#CC79A7","#CC79A7"), bg="transparent", 
       lty=c(2,1), lwd=c(2,3),cex=1,box.lty=0,seg.len = 2,x.intersp=0.5, y.intersp = 0.5)
Fig_S4d <- recordPlot()

pdf("SFig4.pdf",width = 12,height = 12)
plot_grid(Fig_S4a,Fig_S4b,Fig_S4c, Fig_S4d, 
          nrow = 2,labels = c('(a)','(b)','(c)','(d)'))
dev.off()

## Jena
# Model 1: PLs and Ref
jena = filter(Jena, !startsWith(Genotype_ID,"ir")) %>%
  mutate(GeneGroup = ifelse(startsWith(Genotype_ID,"P"),"PL","Ref"))

# Decide main effects on Wavelength 1354, 731, 654, 441
best_model_jena(1354) # M2 < M1
best_model_jena(731)  # M1 < M2
best_model_jena(654) # M1 < M2
best_model_jena(441) # M1 < M2
# We select M1. Genotypic effect would be same as one-way ANOVA

p = vector()
for (i in 1:2101){
  res.aov <- aov(jena[,i+50]~jena$GeneGroup)
  p[i] = summary(res.aov)[[1]][["Pr(>F)"]][1]
}

# check the assumption: residual analysis
plot(res.aov) #OK

# Report model at 500nm
report(aov(jena[,151]~jena$GeneGroup))

p_o = t(as.data.frame(p))
colnames(p_o) = 400:2500
p_a = t(as.data.frame(p.adjust(p, method = "BY")))
colnames(p_a) = 400:2500
sigwavs(p_a) #400-427,502-641,694-1142

plot_quantile(as_spectra(filter(jena,GeneGroup == "Ref")[,51:2151]), col = RGB("#000000"), main="Genotypic effects in Glasshouse",ylab = "Reflectance | P-value",
     ylim=c(0,1.01),total_prob = 1, border=FALSE,xaxs = "i",xlab = "Wavelength (nm)")
plot_quantile(as_spectra(filter(jena,GeneGroup == "PL")[,51:2151]), col=RGB("#E69F00"), total_prob = 1, border=FALSE,add = TRUE)
plot(as_spectra(p_o), col="#CC79A7",lty=2,lwd=2,add = TRUE)
plot(as_spectra(p_a), col="#CC79A7",lwd=3,add = TRUE)
plot_regions(as_spectra(p_o), col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add=TRUE)
abline(v=c(1000,1800), lty = 3)
abline(h=0.05, lty = 3, col = "#CC79A7")
legend(1690,0.4,legend=c("PLs", "Ref","P-value","Adj. p-value"),
       col=c(NA,NA,"#CC79A7","#CC79A7"),fill=c(RGB("#000000"),RGB("#E69F00"),NA,NA), bg="transparent", border="transparent",
       lty=c(NA,NA,2,1),lwd=c(NA,NA,2,3),cex=1,box.lty=0,seg.len = 2,x.intersp=0.5, y.intersp = 0.5)
Fig_6e <- recordPlot()

# Comparison: Ref vs PL in Jena
#Bootstrapping
n = 100
perm = sapply(1:n, function(i){
  boot_id = sample(nrow(Jena_PL),nrow(Jena_Ref),replace=FALSE)
  boot_Jena_PL = Jena_PL[boot_id,51:2151] 
  t(as.data.frame(apply(boot_Jena_PL,2,sd) / apply(boot_Jena_PL,2,mean)))
})

CV_boot_Jena_PL = t(perm)
colnames(CV_boot_Jena_PL) = 400:2500

plot_quantile(as_spectra(CV_boot_Jena_PL), total_prob = 1, border = FALSE, 
              main="Bootstrapped CV, PLs vs Ref, Glasshouse",ylab = "CV",xlab="Wavelength (nm)",
              col = RGB("#E69F00"),xaxs = "i")
plot(mean(as_spectra(CV_boot_Jena_PL)), col="#E69F00",lty = 2,lwd=2,add = TRUE)
plot(as_spectra(CV_Jena_Ref), col="#000000",lwd=3, add = TRUE)
abline(v=c(1000,1800), lty = 3)
plot_regions(as_spectra(CV_Jena_Ref),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add = TRUE)
legend(1750,0.05,legend=c("PLs", "Mean PLs","Ref"),
       col=c(NA,"#E69F00","#000000"),fill=c(RGB("#E69F00"),NA,NA),bg="transparent", border="transparent",
       lty=c(NA,2,1), lwd=c(NA,2,3),cex=1,box.lty=0,seg.len = 2,x.intersp=0.5, y.intersp = 0.5)
Fig_6f <- recordPlot()

pdf("Fig6.pdf",width = 12, height=18)
plot_grid(Fig_6a,Fig_6b,Fig_6c, Fig_6d, Fig_6e, Fig_6f, 
          nrow = 3,labels = c('(a)','(b)','(c)','(d)','(e)','(f)'))
dev.off()

# Model 2: TLs and EVs
jena_t = filter(Jena, Genotype_ID %in% c("pRESC2NC","pSOL3NC") | startsWith(Genotype_ID, "ir"))%>%
  mutate(GeneGroup = ifelse(Genotype_ID %in% c("pRESC2NC","pSOL3NC"),"EV",Genotype_ID))
table(jena_t$Ctime) # am:27, noon:1. Make no sense to include time effect.

p_ACO = vector()
p_AOC = vector()
p_CAD = vector()
p_CHAL = vector()
p_HQT = vector()
p_LOX = vector()
p_RCA = vector()
p_genegroup = vector()

for (i in 1:2101){
  variable = jena_t[,i+50]
  M = lm(variable ~ GeneGroup, data = jena_t)
  p_ACO[i] = summary(M)$coefficients[2,4]
  p_AOC[i] = summary(M)$coefficients[3,4]
  p_CAD[i] = summary(M)$coefficients[4,4]
  p_CHAL[i] = summary(M)$coefficients[5,4]
  p_HQT[i] = summary(M)$coefficients[6,4]
  p_LOX[i] = summary(M)$coefficients[7,4]
  p_RCA[i] = summary(M)$coefficients[8,4]
  A = anova(M)
  p_genegroup[i] = A$`Pr(>F)`[1]
}

# check the assumption: residual analysis
plot(M) #OK

# Report model at 500nm
variable = jena_t[,151]
report(lm(variable ~ GeneGroup, data = jena_t))

ACO = t(as.data.frame(p_ACO))
#ACO = t(as.data.frame(p.adjust(p_ACO, method = "BY"))) # all 1
colnames(ACO) = 400:2500

AOC = t(as.data.frame(p_AOC))
colnames(AOC) = 400:2500
a_AOC = t(as.data.frame(p.adjust(p_AOC, method = "BY"))) # min 0.0012
colnames(a_AOC) = 400:2500
sigwavs(a_AOC) #1135~1889,2034~2381

CAD = t(as.data.frame(p_CAD))
#CAD = t(as.data.frame(p.adjust(p_CAD, method = "BY"))) # all 1
colnames(CAD) = 400:2500

CHAL = t(as.data.frame(p_CHAL)) #min 0.0001726079
#CHAL = t(as.data.frame(p.adjust(p_CHAL, method = "BY"))) # all 1
colnames(CHAL) = 400:2500

HQT = t(as.data.frame(p_HQT)) #min 0.01391814
#HQT = t(as.data.frame(p.adjust(p_HQT, method = "BY"))) # all 1
colnames(HQT) = 400:2500

LOX = t(as.data.frame(p_LOX))
#LOX = t(as.data.frame(p.adjust(p_LOX, method = "BY"))) # all 1
colnames(LOX) = 400:2500

RCA = t(as.data.frame(p_RCA))
a_RCA = t(as.data.frame(p.adjust(p_RCA, method = "BY"))) # min 0.44
colnames(RCA) = 400:2500
colnames(a_RCA) = 400:2500

m6_genogroup_o = t(as.data.frame(p_genegroup))
colnames(m6_genogroup_o) = 400:2500
m6_genogroup_a = t(as.data.frame(p.adjust(p_genegroup, method = "BY")))
colnames(m6_genogroup_a) = 400:2500
sigwavs(m6_genogroup_a) #400-402,1134-1890,2035-2359

plot(as_spectra(AOC), main = "Genotypic effects in Glasshouse", col = "#0072B2", lwd=2,lty=2,
     ylab = "P-value", xlab="Wavelength (nm)",xaxs = "i")
plot(as_spectra(RCA), col = "#E69F00", lwd=2,lty=2,add=TRUE)
#plot(as_spectra(CHAL), col = "#009E73", lwd=1,lty=2,add=TRUE)
#plot(as_spectra(HQT), col = "#F0E442", lwd=2,lty=2,add=TRUE)
plot(as_spectra(a_AOC), col = "#0072B2",lwd=3,add=TRUE)
plot(as_spectra(a_RCA), col = "#E69F00", lwd=3,add=TRUE)
abline(v=c(1000,1800), lty = 3)
abline(h=0.05, lty = 3, col="#CC79A7")
plot_regions(as_spectra(a_AOC),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add=TRUE)
legend(1300,0.35,legend=c("irAOC","Adj. irAOC","irRCA","Adj. irRCA"),
       col=c("#0072B2","#0072B2","#E69F00","#E69F00"), 
       bg="transparent", lty=c(2,1,2,1), lwd=c(2,3,2,3),cex=1,box.lty=0,seg.len = 0.8,x.intersp=0.5,y.intersp = 0.5)
Fig_7c <- recordPlot()

plot(as_spectra(AOC), main = "Genotypic effects in Glasshouse", col = "#0072B2", lwd=2,lty=2,
     ylab = "p-value", xlab="Wavelength (nm)",xaxs = "i",cex.lab=1.4, cex.axis=1.4, cex.main=1.5)
plot(as_spectra(LOX), col = "#D55E00",lwd=2, lty=2,add=TRUE)
plot(as_spectra(RCA), col = "#E69F00", lwd=2,lty=2,add=TRUE)
plot(as_spectra(CHAL), col = "#009E73", lwd=2,lty=2,add=TRUE)
plot(as_spectra(HQT), col = "#F0E442", lwd=2,lty=2,add=TRUE)
plot(as_spectra(CAD), col = "#e31a1c", lwd=2,lty=2,add=TRUE)
plot(as_spectra(ACO), col = "#a6cee3", lwd=2,lty=2,add=TRUE)
plot(as_spectra(a_AOC), col = "#0072B2",lwd=3, add=TRUE)
plot(as_spectra(a_RCA), col = "#E69F00", lwd=3,add=TRUE)
plot_regions(as_spectra(a_AOC),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add=TRUE)
abline(v=c(1000,1800), lty = 3)
abline(h=0.05, lty = 3, col="#CC79A7")
legend(1980,0.5,legend=c("irAOC","Adj.irAOC","irRCA","Adj.irRCA","irLOX3","irCHAL","irHQT","irCAD","irACO"),
       col=c("#0072B2","#0072B2","#E69F00","#E69F00","#D55E00","#009E73","#F0E442","#e31a1c","#a6cee3"), 
       bg="transparent", lty=c(2,1,2,1,2,2,2,2,2), lwd=c(2,3,2,3,2,2,2,2,2),
       cex=1,box.lty=0,seg.len = 0.8,x.intersp=0.5,y.intersp = 0.5)
Fig_S3c <- recordPlot()

#Bootstrapping
n = 100
perm = sapply(1:n, function(i){
  boot_id = sample(nrow(Jena_TL),nrow(Jena_EV),replace=FALSE)
  boot_Jena_TL = Jena_TL[boot_id,51:2151] 
  t(as.data.frame(apply(boot_Jena_TL,2,sd) / apply(boot_Jena_TL,2,mean)))
})

CV_boot_Jena_TL = t(perm)
colnames(CV_boot_Jena_TL) = 400:2500

plot_quantile(as_spectra(CV_boot_Jena_TL), total_prob = 1, border = FALSE, 
              main="Bootstrapped CV, TLs vs. Ref, Glasshouse",ylab = "CV", xlab="Wavelength (nm)",
              col = RGB("#0072B2"),xaxs = "i")
plot(mean(as_spectra(CV_boot_Jena_TL)), col="#0072B2",lty = 2,lwd=2,add = TRUE)
plot(as_spectra(CV_Jena_EV), col="#000000",lwd=3, add = TRUE)
abline(v=c(1000,1800), lty = 3)
plot_regions(as_spectra(CV_Jena_EV),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add = TRUE)
legend(400,0.27,legend=c("TLs", "Mean TLs","EV"),
       col=c(NA,"#0072B2","#000000"), fill=c("#0072B2",NA,NA),bg="transparent",border="transparent", 
       lty=c(NA,2,1), lwd=c(NA,2,3),cex=1,box.lty=0,seg.len = 2,x.intersp=0.5, y.intersp = 0.5)
Fig_7d <- recordPlot()

## Utah
# Model 1: TLs and EV
UT2 = mutate(UT_lin,Plant_ID = paste0("p_",Plant_ID)) 
UT2[UT2$Genotype_ID == "pRESC2NC",]$Genotype_ID = "a_pRESC2NC"
# Decide main effects on Wavelength 1516,1005,536,1801
variable = UT2[,1516-350+1]
M1 <- lme(variable ~ 1 + Ctime + Genotype_ID, data=UT2, random=~1|Plant_ID, method = "REML")
M2 <- lme(variable ~ 1 + Ctime + Genotype_ID, data=UT2, random=~1 + Ctime|Plant_ID, method = "REML")
M3 <- lme(variable ~ 1 + Ctime * Genotype_ID, data=UT2, random=~1|Plant_ID, method = "REML")
M4 <- lme(variable ~ 1 + Ctime * Genotype_ID, data=UT2, random=~1 + Ctime|Plant_ID, method = "REML")
AIC(M1,M2,M3,M4) #M1 < M2 < M3 < M4

variable = UT2[,1005-350+1]
M1 <- lme(variable ~ 1 + Ctime + Genotype_ID, data=UT2, random=~1|Plant_ID, method = "REML")
M2 <- lme(variable ~ 1 + Ctime + Genotype_ID, data=UT2, random=~1 + Ctime|Plant_ID, method = "REML")
M3 <- lme(variable ~ 1 + Ctime * Genotype_ID, data=UT2, random=~1|Plant_ID, method = "REML")
M4 <- lme(variable ~ 1 + Ctime * Genotype_ID, data=UT2, random=~1 + Ctime|Plant_ID, method = "REML")
AIC(M1,M2,M3,M4) #M2 < M1 < M4 < M3

variable = UT2[,536-350+1]
M1 <- lme(variable ~ 1 + Ctime + Genotype_ID, data=UT2, random=~1|Plant_ID, method = "REML")
M2 <- lme(variable ~ 1 + Ctime + Genotype_ID, data=UT2, random=~1 + Ctime|Plant_ID, method = "REML")
M3 <- lme(variable ~ 1 + Ctime * Genotype_ID, data=UT2, random=~1|Plant_ID, method = "REML")
M4 <- lme(variable ~ 1 + Ctime * Genotype_ID, data=UT2, random=~1 + Ctime|Plant_ID, method = "REML")
AIC(M1,M2,M3,M4) #M2 < M1 < M4 < M3

variable = UT2[,1801-350+1]
M1 <- lme(variable ~ 1 + Ctime + Genotype_ID, data=UT2, random=~1|Plant_ID, method = "REML")
M2 <- lme(variable ~ 1 + Ctime + Genotype_ID, data=UT2, random=~1 + Ctime|Plant_ID, method = "REML")
M3 <- lme(variable ~ 1 + Ctime * Genotype_ID, data=UT2, random=~1|Plant_ID, method = "REML")
M4 <- lme(variable ~ 1 + Ctime * Genotype_ID, data=UT2, random=~1 + Ctime|Plant_ID, method = "REML")
AIC(M1,M2,M3,M4) #M1 < M2 < M3 < M4
# We select M1

# check the assumption: residual analysis
qqnorm(M1) #OK
qqnorm(M1, ~ranef(.)) #OK

# Report model at 500nm
variable = UT2[,151]
report(lme(variable ~ 1 + Ctime + Genotype_ID, data=UT2, random=~1|Plant_ID, method = "REML"))

p_time = vector()
amnoon = vector()
ampm = vector()
noonpm = vector()
p_ACO = vector()
p_AOC = vector()
p_CAD = vector()
p_CHAL = vector()
p_HQT = vector()
p_LOX = vector()
p_RCA = vector()

for (i in 1:2101){
  variable = UT2[,i+50]
  M = lme(variable ~ 1 + Ctime + Genotype_ID, data=UT2, random=~1|Plant_ID, method = "REML",control=lmeControl(returnObject=TRUE))
  p_ACO[i] = summary(M)$tTable[,5][4]
  p_AOC[i] = summary(M)$tTable[,5][5]
  p_CAD[i] = summary(M)$tTable[,5][6]
  p_CHAL[i] = summary(M)$tTable[,5][7]
  p_HQT[i] = summary(M)$tTable[,5][8]
  p_LOX[i] = summary(M)$tTable[,5][9]
  p_RCA[i] = summary(M)$tTable[,5][10]
 
  A = anova(M)
  p_time[i] = A$`p-value`[2]
  
  pairwise = emmeans(M, pairwise ~ Ctime)
  pairwise = as.data.frame(pairwise$contrasts %>% summary(infer=TRUE))
  amnoon[i] = pairwise$p.value[1]
  ampm[i] = pairwise$p.value[2]
  noonpm[i] = pairwise$p.value[3]
}

# Genotypic effect
m7_ACO = t(as.data.frame(p_ACO))
#ACO = t(as.data.frame(p.adjust(p_ACO, method = "BY"))) #all 1
colnames(m7_ACO) = 400:2500

m7_AOC = t(as.data.frame(p_AOC))
#AOC = t(as.data.frame(p.adjust(p_AOC, method = "BY"))) # all 1
colnames(m7_AOC) = 400:2500

m7_CAD = t(as.data.frame(p_CAD))
#CAD = t(as.data.frame(p.adjust(p_CAD, method = "BY"))) #all 1
colnames(m7_CAD) = 400:2500

m7_CHAL = t(as.data.frame(p_CHAL))
m7_CHAL_a = t(as.data.frame(p.adjust(p_CHAL, method = "BY"))) 
colnames(m7_CHAL) = 400:2500
colnames(m7_CHAL_a) = 400:2500
sigwavs(m7_CHAL_a) #400-420

m7_HQT = t(as.data.frame(p_HQT))
#HQT = t(as.data.frame(p.adjust(p_HQT, method = "BY"))) # all 1
colnames(m7_HQT) = 400:2500

m7_LOX = t(as.data.frame(p_LOX))
m7_LOX_a = t(as.data.frame(p.adjust(p_LOX, method = "BY"))) # all 1
colnames(m7_LOX) = 400:2500
colnames(m7_LOX_a) = 400:2500

m7_RCA = t(as.data.frame(p_RCA))
#RCA = t(as.data.frame(p.adjust(p_RCA, method = "BY"))) # all 1
colnames(m7_RCA) = 400:2500

plot(as_spectra(m7_CHAL), main = "Genotypic effects in Field_UT", ylim = c(0,1.01), xaxs = "i",
     col = "#009E73", lty=2,lwd=2,ylab = "P-value", xlab="Wavelength (nm)")
plot(as_spectra(m7_LOX), col = "#D55E00", lwd=2,lty=2,add=TRUE)
plot(as_spectra(m7_CHAL_a), col = "#009E73",lwd=3, add=TRUE)
plot(as_spectra(m7_LOX_a), col = "#D55E00", lwd=3,add=TRUE)
plot_regions(as_spectra(m7_CHAL_a),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add=TRUE)
abline(v=c(1000,1800), lty = 3)
abline(h=0.05, lty = 3, col="#CC79A7")
legend(1800,0.3,legend=c("irCHAL","Adj. irCHAL","irLOX3","Adj. irLOX3"),
       col=c("#009E73","#009E73","#D55E00","#D55E00"), 
       bg="transparent", lty=c(2,1,2,1), lwd=c(2,3,2,3),cex=1,box.lty=0,seg.len = 2,x.intersp=0.5,y.intersp = 0.5)
Fig_7a <- recordPlot()

plot(as_spectra(m7_AOC), main = "Genotypic effects in Field_UT", ylim = c(0,1.01), xaxs = "i",
     col = "#0072B2", lty=2,lwd=2,ylab = "P-value", xlab="Wavelength (nm)",cex.lab=1.4, cex.axis=1.4, cex.main=1.5)
plot(as_spectra(m7_LOX), col = "#D55E00",lty=2,lwd=2, add=TRUE)
plot(as_spectra(m7_RCA), col = "#E69F00", lty=2,lwd=2,add=TRUE)
plot(as_spectra(m7_CHAL), col = "#009E73", lty=2,lwd=2,add=TRUE)
plot(as_spectra(m7_HQT), col = "#F0E442", lty=2,lwd=2,add=TRUE)
plot(as_spectra(m7_CAD), col = "#e31a1c", lty=2,lwd=2,add=TRUE)
plot(as_spectra(m7_ACO), col = "#a6cee3", lty=2,lwd=2,add=TRUE)
plot(as_spectra(m7_CHAL_a), col = "#009E73",lwd=3, add=TRUE)
plot(as_spectra(m7_LOX_a), col = "#D55E00",lwd=3, add=TRUE)
plot_regions(as_spectra(m7_AOC),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add=TRUE)
abline(v=c(1000,1800), lty = 3)
abline(h=0.05, lty = 3, col="#CC79A7")
legend(1900,0.38,legend=c("irCHAL","Adj. irCHAL","irLOX3","Adj. irLOX3","irAOC","irRCA","irHQT","irCAD","irACO"),
       col=c("#009E73","#009E73","#D55E00","#D55E00","#0072B2","#E69F00","#F0E442","#e31a1c","#a6cee3"), 
       bg="transparent", lty=c(2,1,2,1,2,2,2,2,2), lwd=c(2,3,2,3,2,2,2,2,2),
       cex=1,box.lty=0,seg.len = 0.8,x.intersp=0.5,y.intersp = 0.5)
Fig_S3b <- recordPlot()

pdf("SFig3.pdf",width = 12,height = 12)
plot_grid(Fig_S3a,NULL,Fig_S3b,Fig_S3c, 
          nrow = 2,labels = c('(a)','','(b)','(c)'))
dev.off()

#Bootstrapping
n = 100
perm = sapply(1:n, function(i){
  boot_id = sample(nrow(UT_TL),nrow(UT_EV),replace=FALSE)
  boot_UT_TL = UT_TL[boot_id,51:2151] 
  t(as.data.frame(apply(boot_UT_TL,2,sd) / apply(boot_UT_TL,2,mean)))
})

CV_boot_UT_TL = t(perm)
colnames(CV_boot_UT_TL) = 400:2500

plot_quantile(as_spectra(CV_boot_UT_TL), total_prob = 1, border = FALSE, 
              main="Bootstrapped CV, TLs vs Ref, Field_UT",ylab = "Coefficient of variance",xlab="Wavelength (nm)",
              col = RGB("#0072B2"),xaxs = "i")
plot(mean(as_spectra(CV_boot_UT_TL)), col="#0072B2",lty = 2,lwd=2,add = TRUE)
plot(as_spectra(CV_UT_EV), col="#000000",lwd=3, add = TRUE)
abline(v=c(1000,1800), lty = 3)
plot_regions(as_spectra(CV_UT_EV),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add = TRUE)
legend(400,0.4,legend=c("TLs", "Mean TLs","Ref"),
       col=c(NA,"#0072B2","#000000"), fill=c("#0072B2",NA,NA),bg="transparent",border="transparent", 
       lty=c(NA,2,1), lwd=c(NA,2,3),cex=1,box.lty=0,seg.len = 2,x.intersp=0.5, y.intersp = 0.5)
Fig_7b <- recordPlot()

pdf("Fig7.pdf",width = 12,height = 12)
plot_grid(Fig_7a,Fig_7b,Fig_7c,Fig_7d,
          nrow=2,labels = c('(a)','(b)','(c)','(d)'))
dev.off()


# Time effect
m7_p_time_o = t(as.data.frame(p_time))
m7_p_time_a = t(as.data.frame(p.adjust(p_time, method = "BY")))
colnames(m7_p_time_o) = 400:2500 
colnames(m7_p_time_a) = 400:2500 
sigwavs(m7_p_time_a) #424~611,696-1364,1511-1788

plot_quantile(as_spectra(filter(UT2,Ctime=="am")[,51:2151]),ylim = c(0,1.01), total_prob = 1, border=FALSE,col = RGB("#F0E442"), main = "Effects of time in Field_UT", 
              ylab = "Reflectance | P-value", xlab="Wavelength (nm)",xaxs = "i")
plot_quantile(as_spectra(filter(UT2,Ctime=="noon")[,51:2151]),col = RGB("#E69F00"), total_prob = 1, border=FALSE,add=TRUE)
plot_quantile(as_spectra(filter(UT2,Ctime=="pm")[,51:2151]),col = RGB("#0072B2"), total_prob = 1, border=FALSE,add=TRUE)
plot(as_spectra(m7_p_time_o),col = "#CC79A7", lwd=2, lty=2, add=TRUE)
plot(as_spectra(m7_p_time_a), col = "#CC79A7", lwd=3, add=TRUE)
abline(v=c(1000,1800), lty = 3)
abline(h=0.05, lty = 3, col="#CC79A7")
plot_regions(as_spectra(m7_p_time_o),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add=TRUE)
legend(400,1,legend=c("Am", "Noon","Pm","P-value","Adj. p-value"),
       fill=c(RGB("#F0E442"),RGB("#E69F00"),RGB("#0072B2"),NA,NA), col=c(NA,NA,NA,"#CC79A7","#CC79A7"),bg="transparent", border="transparent",
       lty=c(NA,NA,NA,2,1), lwd=c(NA,NA,NA,2,3),cex=1,box.lty=0,seg.len = 2,x.intersp=0.5, y.intersp = 0.5)
Fig_8a <-  recordPlot()

m7_amnoon_o = t(as.data.frame(amnoon))
colnames(m7_amnoon_o) = 400:2500
m7_amnoon_a = t(as.data.frame(p.adjust(amnoon, method = "BY")))
colnames(m7_amnoon_a) = 400:2500

m7_ampm_o = t(as.data.frame(ampm))
colnames(m7_ampm_o) = 400:2500
m7_ampm_a = t(as.data.frame(p.adjust(ampm, method = "BY")))
colnames(m7_ampm_a) = 400:2500
sigwavs(m7_ampm_a) 

m7_noonpm_o = t(as.data.frame(noonpm))
colnames(m7_noonpm_o) = 400:2500
m7_noonpm_a = t(as.data.frame(p.adjust(noonpm, method = "BY")))
colnames(m7_noonpm_a) = 400:2500
sigwavs(m7_noonpm_a) # 425-613,696-1114,1163-1359,1511-1785

plot(as_spectra(m7_amnoon_o), col="#F0E442",main="Post-hoc,Effect of time in Field_UT",xlab="Wavelength (nm)", ylab = "P-value",
     xaxs = "i",lty=2,lwd=2,ylim=c(0,1.01))
plot(as_spectra(m7_ampm_o), col="#0072B2", lty=2,lwd=2,add = TRUE)
plot(as_spectra(m7_noonpm_o), col="#E69F00",lty=2, lwd=2,add = TRUE)
plot(as_spectra(m7_amnoon_a), col="#F0E442", lwd=3,add = TRUE)
plot(as_spectra(m7_ampm_a), col="#0072B2", lwd=3,add = TRUE)
plot(as_spectra(m7_noonpm_a), col="#E69F00",lwd=3, add = TRUE)
abline(v=c(1000,1800), lty = 3)
abline(h=0.05, lty = 3, col = "#CC79A7")
plot_regions(as_spectra(m7_ampm_o),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add = TRUE)
legend(1800,0.2,legend=c("Am vs Noon","AM vs Pm","Noon vs Pm"),
       col=c("#F0E442","#0072B2","#E69F00"), bg="transparent", 
       lty=1, lwd=3,cex=1,box.lty=0,seg.len = 2,x.intersp=0.5, y.intersp = 0.5)
Fig_8b <- recordPlot()

pdf("Fig8.pdf",height = 12,width = 12)
plot_grid(Fig_8a,Fig_8b,Fig_8c, Fig_8d, 
          nrow = 2,labels = c('(a)','(b)','(c)','(d)'))
dev.off()

pdf("Fig9.pdf",width = 12,height = 12)
plot_grid(Fig_9a,Fig_9b,Fig_9c, Fig_9d, 
          nrow = 2,labels = c('(a)','(b)','(c)','(d)'))
dev.off()

# Model 2: UT_dev
# UT_dev leaf position
utdev = mutate(UT_dev,leaf = ifelse(Leaf_ID == 2,"a2",paste0("b",Leaf_ID)),
               Plant_ID = paste0("p",Plant_ID))

# Decide main effects on Wavelength 1398,689,718,711
variable = utdev[,1398-350+1]
M1 <- lme(variable ~ 1 + Ctime + leaf, data=utdev, random=~1|Plant_ID, method = "REML")
M2 <- lme(variable ~ 1 + Ctime * leaf, data=utdev, random=~1|Plant_ID, method = "REML")
AIC(M1,M2) #M1 < M2

variable = utdev[,689-350+1]
M1 <- lme(variable ~ 1 + Ctime + leaf, data=utdev, random=~1|Plant_ID, method = "REML")
M2 <- lme(variable ~ 1 + Ctime * leaf, data=utdev, random=~1|Plant_ID, method = "REML")
AIC(M1,M2) #M1 < M2 

variable = utdev[,718-350+1]
M1 <- lme(variable ~ 1 + Ctime + leaf, data=utdev, random=~1|Plant_ID, method = "REML")
M2 <- lme(variable ~ 1 + Ctime * leaf, data=utdev, random=~1|Plant_ID, method = "REML")
AIC(M1,M2) #M1 < M2

variable = utdev[,711-350+1]
M1 <- lme(variable ~ 1 + Ctime + leaf, data=utdev, random=~1|Plant_ID, method = "REML")
M2 <- lme(variable ~ 1 + Ctime * leaf, data=utdev, random=~1|Plant_ID, method = "REML")
AIC(M1,M2) #M1 < M2
# We select M1

# check the assumption: residual analysis
qqnorm(M1) #OK
qqnorm(M1, ~ranef(.)) #OK

# Report model at 500nm
variable = utdev[,151]
report(lme(variable ~ 1 + Ctime + leaf, data=utdev, random=~1|Plant_ID, method = "REML"))

p_noon = vector()
avo_leaf = vector()
l12 = vector()
l23 = vector()
l13 = vector()

for (i in 1:2101){
  variable = utdev[,i+50]
  M = lme(variable ~ 1 + Ctime + leaf, data=utdev, random=~1|Plant_ID, method = "REML")
  p_noon[i] = summary(M)$tTable[,5][2]

  A = anova(M)
  avo_leaf[i] = A$`p-value`[3]
  
  pairwise = emmeans(M, pairwise ~ leaf)
  pairwise = as.data.frame(pairwise$contrasts %>% summary(infer=TRUE))
  l12[i] = pairwise$p.value[1]
  l23[i] = pairwise$p.value[2]
  l13[i] = pairwise$p.value[3]
}

# Effect of leaf position
l12_o = t(as.data.frame(l12))
l12_a = t(as.data.frame(p.adjust(l12, method = "BY")))
colnames(l12_o) = 400:2500
colnames(l12_a) = 400:2500

l23_o = t(as.data.frame(l23))
l23_a = t(as.data.frame(p.adjust(l23, method = "BY")))
colnames(l23_o) = 400:2500
colnames(l23_a) = 400:2500

l13_o = t(as.data.frame(l13))
l13_a = t(as.data.frame(p.adjust(l13, method = "BY")))
colnames(l13_o) = 400:2500
colnames(l13_a) = 400:2500

avo_leaf_o = t(as.data.frame(avo_leaf))
avo_leaf_a = t(as.data.frame(p.adjust(avo_leaf, method = "BY")))
colnames(avo_leaf_o) = 400:2500
colnames(avo_leaf_a) = 400:2500 
sigwavs(avo_leaf_a) 
as.numeric(rownames(t(avo_leaf_a))[t(avo_leaf_a)<1])

plot_quantile(as_spectra(filter(utdev,leaf=="b1")[,51:2151]),ylim = c(0,1.01), total_prob = 1, border=FALSE,col = RGB("#F0E442"),main = "Effects of leaf position in Field_UT", 
              xaxs = "i",ylab = "Reflectance | P-value", xlab="Wavelength (nm)")
plot_quantile(as_spectra(filter(utdev,leaf=="a2")[,51:2151]),col = RGB("#E69F00"), total_prob = 1, border=FALSE,add=TRUE)
plot_quantile(as_spectra(filter(utdev,leaf=="b3")[,51:2151]),col = RGB("#0072B2"), total_prob = 1, border=FALSE,add=TRUE)
plot(as_spectra(avo_leaf_o),col = "#CC79A7", lwd=2, lty=2, add=TRUE)
plot(as_spectra(avo_leaf_a), col = "#CC79A7", lwd=3, add=TRUE)
abline(v=c(1000,1800), lty = 3)
abline(h=0.05, lty = 3, col="#CC79A7")
plot_regions(as_spectra(avo_leaf_o),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add=TRUE)
legend(1800,0.5,legend=c("Leaf 1", "Leaf 2","Leaf 3","P-value","Adj. p-value"),
       fill=c(RGB("#F0E442"),RGB("#E69F00"),RGB("#0072B2"),NA,NA), col=c(NA,NA,NA,"#CC79A7","#CC79A7"),bg="transparent", border="transparent",
       lty=c(NA,NA,NA,2,1), lwd=c(NA,NA,NA,2,3),cex=1,box.lty=0,seg.len = 2,x.intersp=0.5, y.intersp = 0.5)
Fig_9c <-  recordPlot()

plot(as_spectra(l12_o), col="#F0E442",main="Effect of leaf position in Field_UT, Post-hoc",xlab="Wavelength (nm)", ylab = "P-value",
     xaxs = "i",lty=2,lwd=2,ylim=c(0,1.01))
plot(as_spectra(l13_o), col="#0072B2", lty=2,lwd=2,add = TRUE)
plot(as_spectra(l23_o), col="#E69F00",lty=2, lwd=2,add = TRUE)
plot(as_spectra(l12_a), col="#F0E442", lwd=3,add = TRUE)
plot(as_spectra(l13_a), col="#0072B2", lwd=3,add = TRUE)
plot(as_spectra(l23_a), col="#E69F00",lwd=3, add = TRUE)
abline(v=c(1000,1800), lty = 3)
abline(h=0.05, lty = 3, col = "#CC79A7")
plot_regions(as_spectra(l12_a),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add = TRUE)
legend(1850,0.4,legend=c("L1 vs L2","L1 vs L3","L2 vs L3"),
       col=c("#F0E442","#0072B2","#E69F00"), bg="transparent", 
       lty=1, lwd=3,cex=1,box.lty=0,seg.len = 2,x.intersp=0.5, y.intersp = 0.5)
Fig_9d <- recordPlot()

# Effect of time
m8_p_noon_o = t(as.data.frame(p_noon))
m8_p_noon_a = t(as.data.frame(p.adjust(p_noon, method = "BY")))
colnames(m8_p_noon_o) = 400:2500
colnames(m8_p_noon_a) = 400:2500
sigwavs(m8_p_noon_a) #None

plot(as_spectra(m8_p_noon_o), main = "Effects of time in Field_UT", ylim = c(0,1.01),col = "#CC79A7", lty=2, lwd=2, ylab = "P-value", xlab="Wavelength (nm)",xaxs = "i")
plot(as_spectra(m8_p_noon_a), lwd=3, col = "#CC79A7", add=TRUE)
abline(v=c(1000,1800), lty = 3)
abline(h=0.05, lty = 3, col="#CC79A7")
plot_regions(as_spectra(m8_p_noon_a),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add=TRUE)
legend(1800,0.3,legend=c("P-value","Adj. p-value"),
       col=c("#CC79A7","#CC79A7"), bg="transparent", 
       lty=c(2,1), lwd=c(2,3),cex=1,box.lty=0,seg.len = 2,x.intersp=0.5, y.intersp = 0.5)
Fig_S4c <- recordPlot()

pdf("SFig4.pdf",width = 12,height = 12)
plot_grid(Fig_S4a,Fig_S4b,Fig_S4c,Fig_S4d,
          nrow=2,labels = c('(a)','(b)','(c)','(d)'))
dev.off()


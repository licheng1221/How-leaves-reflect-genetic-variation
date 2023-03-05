library(MASS)
library(tidyverse)
library(spectrolab)
library(lme4)
library(nlme)
library(emmeans)
library(report)

# Load customized functions
source("~/Desktop/Nicotiana/Final/Code/Functions.R")
options(max.print = 3000)

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

# Report model at 500nm
report(aov(All_UT[,101]~All_UT$Site))

p_o = t(as.data.frame(p))
colnames(p_o) = 400:2500
p_a = t(as.data.frame(p.adjust(p, method = "BY")))
colnames(p_a) = 400:2500

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

par(mar=c(5.5, 5.5,4, 6) + 0.1)
plot(as_spectra(UT_UT[,51:2101]), col="#FC4E07", main="Reference (UT-WT) across experiments",xlab="Wavelength", ylab = "Reflectance",
     cex.lab=1.4, cex.axis=1.4, cex.main=1.5, ylim=c(0,1.01), xaxs = "i")
plot(as_spectra(AZ_UT[,51:2101]), col="#00AFBB", add = TRUE)
plot(as_spectra(Jena_UT[,51:2101]), col="#E7B800", add = TRUE)
plot(as_spectra(CV_UT_UT), col="#FC4E07",lty=2, lwd=2, add = TRUE)
plot(as_spectra(CV_Jena_UT), col="#E7B800", lty=2, lwd=2,add = TRUE)
plot(as_spectra(CV_AZ_UT), col="#00AFBB", lty=2, lwd=2,add = TRUE)
plot(as_spectra(p_o),col = "#a50000",lty=3, add=TRUE)
plot(as_spectra(p_a),col = "#a50000",add=TRUE)
abline(v=c(1000,1800), lty = 3)
abline(h=0.05, lty = 3, col = "#a50000")
axis(4, ylim=c(0,1.01), cex.axis=1.4)
plot_regions(as_spectra(CV_UT_UT), col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add=TRUE)
mtext("Coefficient of Variance/P-value",side=4,line=3,cex = 1.5) 
legend(1850,1,legend=c("CV Field_AZ", "CV Field_UT","CV Glasshouse","Adjusted p-value"),
       col=c("#00AFBB","#FC4E07","#E7B800","#a50000"), bg="transparent", 
       lty=c(2,2,2,1), cex=0.8,box.lty=0,seg.len = 1,x.intersp=0.5, y.intersp=1.2)

par(mar=c(5.5, 5.5,4,4))
plot(as_spectra(azut_o), col="grey60", main="Post hoc pairwise caomparison",xlab="Wavelength", ylab = "p-value",
     xaxs = "i", cex.lab=1.4, cex.axis=1.4, lty=2,cex.main=1.5, ylim=c(0,1.01))
plot(as_spectra(azjena_o), col="#00AFBB", lty=2,add = TRUE)
plot(as_spectra(utjena_o), col="#FC4E07", lty=2,add = TRUE)
plot(as_spectra(azut_a), col="grey60", lwd=2,add = TRUE)
plot(as_spectra(azjena_a), col="#00AFBB", lwd=2,add = TRUE)
plot(as_spectra(utjena_a), col="#FC4E07",lwd=2, add = TRUE)
abline(v=c(1000,1800), lty = 3)
abline(h=0.05, lty = 3, col = "#a50000")
plot_regions(as_spectra(p_o),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add = TRUE)
legend(1880,1,legend=c("Field_AZ vs Glasshouse","Field_UT vs Glasshouse","Field_AZ vs Field_UT"),
       col=c("#00AFBB","#FC4E07","grey60"), bg="transparent", 
       lty=c(1,1,1), cex=0.7,box.lty=0,seg.len = 1,x.intersp=0.5,y.intersp = 1.2)

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

p_o = t(as.data.frame(p))
colnames(p_o) = 400:2500
p_a = t(as.data.frame(p.adjust(p, method = "BY")))
colnames(p_a) = 400:2500

par(mar=c(5.5, 5.5,4, 6) + 0.1)
plot(as_spectra(Jena_UT[,1:2101]), col = "grey60", main="Reference samples in Glasshouse",ylab = "Reflectance",
     xaxs = "i",cex.lab=1.4, cex.axis=1.4, cex.main=1.5, ylim=c(0,1.01),xaxs = "i",xlab = "Wavelength")
plot(as_spectra(Jena_EV1[,1:2101]), col="#B8860B",alpha = 0.1, add = TRUE)
plot(as_spectra(Jena_EV2[,1:2101]), col="#F0E68C", add = TRUE)
plot(as_spectra(CV_Jena_UT), col="grey60", lty=2, lwd=2,add = TRUE)
plot(as_spectra(CV_Jena_EV1), col="#B8860B",lty=2, lwd=2, add = TRUE)
plot(as_spectra(CV_Jena_EV2), col="#F0E68C",lty=2, lwd=2, add = TRUE)
plot(as_spectra(p_o), col="#a50000",lty=2,add = TRUE)
plot(as_spectra(p_a), col="#a50000",add = TRUE)
plot_regions(as_spectra(Jena_UT[,1:2101]), col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add=TRUE)
abline(v=c(1000,1800), lty = 3)
abline(h=0.05, lty = 3, col = "#a50000")
axis(4, ylim=c(0,1.01), cex.axis=1.4)
mtext("Coefficient of Variance/ p-value",side=4,line=3,cex = 1.5) 
legend(1900,1,legend=c("UT-WT", "EV1","EV2","P-value"),
       col=c("grey60","#B8860B","#F0E68C","#a50000"), bg="transparent", 
       lty=c(1,1,1,1), cex=0.8,box.lty=0,seg.len = 1,x.intersp=0.5, y.intersp=1.2)

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
# Adding interaction won't help. So we stay with M3

# check the assumption: residual analysis
plot(lm(az[,1877-350+1] ~ GeneGroup + Ctime + Batch + Leaf_Num, data = az)) #OK

# Run linear regression across all wavelengths
p_genegroup = vector()
plril = vector()
plwt = vector()
rilwt = vector()
p_time = vector()
p_noon = vector()
p_pm = vector()
p_batch = vector()
p_leafnum = vector()

for (i in 1:2101){
  variable = az[,i+50]
  M = lm(variable ~ GeneGroup + Ctime + Batch + Leaf_Num, data = az) 
  p_noon[i] = summary(M)$coefficients[4,4]
  p_pm[i] = summary(M)$coefficients[5,4]
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
}

# Report model at 500nm
variable = az[,151]
report(lm(variable ~ GeneGroup + Ctime + Batch + Leaf_Num, data = az))

# Genotypic effects
genegroup_o = t(as.data.frame(p_genegroup))
genegroup_a = t(as.data.frame(p.adjust(p_genegroup, method = "BY")))
colnames(genegroup_o) = 400:2500 
sigwavs(genegroup_o) #400~1900
colnames(genegroup_a) = 400:2500 
sigwavs(genegroup_a) #409~760, 1001~1017,1128-1800

plril_o = t(as.data.frame(plril))
colnames(plril_o) = 400:2500
plril_a = t(as.data.frame(p.adjust(plril, method = "BY")))
colnames(plril_a) = 400:2500

plwt_o = t(as.data.frame(plwt))
colnames(plwt_o) = 400:2500
plwt_a = t(as.data.frame(p.adjust(plwt, method = "BY")))
colnames(plwt_a) = 400:2500
sigwavs(plwt_a) #466-656

rilwt_o = t(as.data.frame(rilwt))
colnames(rilwt_o) = 400:2500
rilwt_a = t(as.data.frame(p.adjust(rilwt, method = "BY")))
colnames(rilwt_a) = 400:2500
sigwavs(rilwt_a) #413-1528

par(mar=c(5.5, 5.5,4, 6) + 0.1)
plot(as_spectra(filter(az,GeneGroup == "RIL")[,51:2151]), col = "#a6cee3", main="Genotypic effects in Field_AZ",ylab = "Reflectance",
     cex.lab=1.4, cex.axis=1.4, cex.main=1.5, ylim=c(0,1.01),xaxs = "i",xlab = "Wavelength")
plot(as_spectra(filter(az,GeneGroup == "aUT-WT")[,51:2151]), col="grey60", add = TRUE)
plot(as_spectra(filter(az,GeneGroup == "PL")[,51:2151]), col="#1f78b4",add = TRUE)
plot(as_spectra(CV_AZ_RIL), col="#a6cee3", lty=2, lwd=2,add = TRUE)
plot(as_spectra(CV_AZ_PL), col="#1f78b4",lty=2, lwd=2, add = TRUE)
plot(as_spectra(CV_AZ_UT), col="grey60",lty=2, lwd=2, add = TRUE)
plot(as_spectra(genegroup_o), col="#a50000",lty=2,lwd=2,add = TRUE)
plot(as_spectra(genegroup_a), col="#a50000",lwd=2,add = TRUE)
plot_regions(as_spectra(CV_AZ_RIL), col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add=TRUE)
abline(v=c(1000,1800), lty = 3)
abline(h=0.05, lty = 3, col = "#a50000")
axis(4, ylim=c(0,1.01), cex.axis=1.4)
mtext("Coefficient of Variance/ p-value",side=4,line=3,cex = 1.5) 
legend(400,1.05,legend=c("PLs", "RILs","UT-WT","p-value"),
       col=c("#1f78b4","#a6cee3","grey60","#a50000"), bg="transparent", 
       lty=c(1,1,1,1), cex=0.8,box.lty=0,seg.len = 1,x.intersp=0.5, y.intersp=1.2)

plot(as_spectra(plril_o), col="grey60",main="Post hoc pairwise comparision",xlab="Wavelength", ylab = "P-value",
     cex.lab=1.4, cex.axis=1.4, xaxs = "i",lty=2,lwd=2,cex.main=1.5, ylim=c(0,1.01))
plot(as_spectra(plwt_o), col="#1f78b4", lty=2,lwd=2,add = TRUE)
plot(as_spectra(rilwt_o), col="#a6cee3",lty=2, lwd=2,add = TRUE)
plot(as_spectra(plril_a), col="grey60", lwd=2,add = TRUE)
plot(as_spectra(plwt_a), col="#1f78b4", lwd=2,add = TRUE)
plot(as_spectra(rilwt_a), col="#a6cee3",lwd=2, add = TRUE)
abline(v=c(1000,1800), lty = 3)
abline(h=0.05, lty = 3, col = "#a50000")
plot_regions(as_spectra(CV_AZ_UT),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add = TRUE)
legend(1900,0.1,legend=c("PLs vs UT-WT","RILs vs UT-WT","PLs vs RILs"),
       col=c("#1f78b4","#a6cee3","grey60"), bg="transparent", 
       lty=c(1,1,1), cex=0.8,box.lty=0,seg.len = 1,x.intersp=0.5, y.intersp=1.2)

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

par(mar=c(5.5, 5.5,4, 4))
plot_quantile(as_spectra(CV_boot_AZ_RIL), total_prob = 1, border = FALSE, 
              main="Bootstrapped CV of RILs vs reference (UT-WT) in Field_AZ",ylab = "Coefficient of variance",xlab="Wavelength",
              col = "#a6cee3",cex.lab=1.4, cex.axis=1.4, cex.main=1.5, xaxs = "i")
plot(mean(as_spectra(CV_boot_AZ_RIL)), col="#1f78b4",lty = 2,add = TRUE)
plot(as_spectra(CV_AZ_UT), col="grey60",lwd=2, add = TRUE)
abline(v=c(1000,1800), lty = 3)
plot_regions(as_spectra(CV_AZ_UT),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add = TRUE)
legend(1950,0.215,legend=c("CV RILs", "Mean CV RILs","CV UT-WT"),
       col=c("#a6cee3","#1f78b4","grey60"), bg="transparent", 
       lty=c(1,2,1,2), cex=0.8,box.lty=0,seg.len = 1,x.intersp=0.5, y.intersp=1.2)

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
              main="CV of PLs vs reference (UT-WT) in Field_AZ",ylab = "Coefficient of variance",xlab="Wavelength",
              col = "grey60",cex.lab=1.4, cex.axis=1.4, cex.main=1.5,xaxs = "i")
plot(mean(as_spectra(CV_boot_AZ_UT)), lty = 2,add = TRUE)
plot(as_spectra(CV_AZ_PL), col="#1f78b4",lwd=2, add = TRUE)
abline(v=c(1000,1800), lty = 3)
plot_regions(as_spectra(CV_AZ_PL),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add = TRUE)
legend(1930,0.06,legend=c("CV PLs","CV UT-WT","Mean CV UT-WT"),
       col=c("#1f78b4","grey60","black"), bg="transparent", 
       lty=c(1,1,1), cex=0.8,box.lty=0,seg.len = 1,x.intersp=0.5, y.intersp=1.2)

# Time effect
noon_o = t(as.data.frame(p_noon))
noon_a = t(as.data.frame(p.adjust(p_noon, method = "BY")))
colnames(noon_o) = 400:2500
colnames(noon_a) = 400:2500

pm_o = t(as.data.frame(p_pm))
pm_a = t(as.data.frame(p.adjust(p_pm, method = "BY")))
colnames(pm_o) = 400:2500
colnames(pm_a) = 400:2500

plot_quantile(as_spectra(filter(az,Ctime=="am")[,51:2151]),ylim = c(0,1.01), total_prob = 1, border=FALSE,col = rgb(0.992, 0.749, 0.435,0.5),main = "Effects of time in Field_AZ", 
              ylab = "Reflectance", xlab="Wavelength",xaxs = "i",cex.lab=1.4, cex.axis=1.4, cex.main=1.5)
plot_quantile(as_spectra(filter(az,Ctime=="noon")[,51:2151]),col = rgb(1, 0.498, 0,0.5), total_prob = 1, border=FALSE,add=TRUE)
plot_quantile(as_spectra(filter(az,Ctime=="pm")[,51:2151]),col = rgb(0.122, 0.471, 0.706,0.5), total_prob = 1, border=FALSE,add=TRUE)
plot(as_spectra(noon_o),col = "#ff7f00", lwd=2, lty=2, add=TRUE)
plot(as_spectra(noon_a), col = "#ff7f00", lwd=2, add=TRUE)
plot(as_spectra(pm_o), col = "#1f78b4", lwd = 2,lty=2, add=TRUE)
plot(as_spectra(pm_a), col = "#1f78b4", lwd=2, add=TRUE)
abline(v=c(1000,1800), lty = 3)
abline(h=0.05, lty = 3, col="red")
plot_regions(as_spectra(noon_o),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add=TRUE)
axis(side = 4, cex.axis=1.5)
mtext("p-value",side=4,line=3,cex = 1.5) 
legend(2050,1,legend=c("batch", "Time","Leaf number"),
       col=c("#ff7f00","#1f78b4","#33a02c"), bg="transparent", 
       lty=c(1,1,1), cex=0.8,box.lty=0,seg.len = 1,x.intersp=0.5, y.intersp=1.2)

# Effect of leaf numbers
leafnum_o = t(as.data.frame(p_leafnum))
leafnum_a = t(as.data.frame(p.adjust(p_leafnum, method = "BY")))
colnames(leafnum_o) = 400:2500
colnames(leafnum_a) = 400:2500 
sigwavs(leafnum_a) # 508~641, 693~707,731~1149,1400~1539,2344~2500

plot_quantile(as_spectra(filter(az,Leaf_Num=="n_1")[,51:2151]),ylim = c(0,1.01), total_prob = 1, border=FALSE,col = rgb(0.992, 0.749, 0.435,0.5),main = "Effects of 1 vs. 2 leaves in Field_AZ", 
              ylab = "Reflectance", xlab="Wavelength",xaxs = "i",cex.lab=1.4, cex.axis=1.4, cex.main=1.5)
plot_quantile(as_spectra(filter(az,Leaf_Num=="n_2")[,51:2151]),col = rgb(0.122, 0.471, 0.706,0.5), total_prob = 1, border=FALSE,add=TRUE)
plot(as_spectra(leafnum_o),col = "#a50000", lwd=2, lty=2, add=TRUE)
plot(as_spectra(leafnum_a), col = "#a50000", lwd=2, add=TRUE)
abline(v=c(1000,1800), lty = 3)
abline(h=0.05, lty = 3, col="red")
plot_regions(as_spectra(leafnum_o),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add=TRUE)
axis(side = 4, cex.axis=1.5)
mtext("p-value",side=4,line=3,cex = 1.5) 

# Effect of batches
batch_o = t(as.data.frame(p_batch))
batch_a = t(as.data.frame(p.adjust(p_batch, method = "BY")))
colnames(batch_o) = 400:2500
colnames(batch_a) = 400:2500 # >0.05: 1898~2001

plot(as_spectra(batch_o), main = "Effects of blocks in Field_AZ", ylim = c(0,0.1),col = "#a50000", lty=3, lwd=2, ylab = "P-value", xlab="Wavelength",
     xaxs = "i",cex.lab=1.4, cex.axis=1.4, cex.main=1.5)
plot(as_spectra(batch_a), lwd=2, col = "#a50000", add=TRUE)
abline(v=c(1000,1800), lty = 3)
abline(h=0.05, lty = 3, col="red")
plot_regions(as_spectra(batch_o),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add=TRUE)

## Jena
# Model 1: PLs and Ref
jena = filter(Jena, !startsWith(Genotype_ID,"ir")) %>%
  mutate(GeneGroup = ifelse(startsWith(Genotype_ID,"P"),"PL","Ref"))

# Decide main effects on Wavelength 1353, 731, 650, 534
best_model_jena(1353) # M2 < M1
best_model_jena(731)  # M1 < M2
best_model_jena(650) # M1 < M2
best_model_jena(534) # M1 < M2
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

par(mar=c(5.5, 5.5,4, 6) + 0.1)
plot(as_spectra(filter(jena,GeneGroup == "Ref")[,51:2151]), col = "grey60", main="PLs, Ref (UT-WT, EVs) in Glasshouse",ylab = "Reflectance",
     cex.lab=1.4, cex.axis=1.4, cex.main=1.5, ylim=c(0,1.01),xaxs = "i",xlab = "Wavelength")
plot(as_spectra(filter(jena,GeneGroup == "PL")[,51:2151]), col="#1f78b4", add = TRUE)
plot(as_spectra(CV_Jena_PL), col="#1f78b4",lty=2, lwd=2, add = TRUE)
plot(as_spectra(CV_Jena_Ref), col="grey60",lty=2, lwd=2, add = TRUE)
plot(as_spectra(p_o), col="#a50000",lty=2,add = TRUE)
plot(as_spectra(p_a), col="#a50000",add = TRUE)
plot_regions(as_spectra(CV_Jena_PL), col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add=TRUE)
abline(v=c(1000,1800), lty = 3)
abline(h=0.05, lty = 3, col = "#a50000")
axis(4, ylim=c(0,1.01), cex.axis=1.4)
mtext("Coefficient of Variance/P-value",side=4,line=3,cex = 1.5) 
legend(400,1.05,legend=c("PLs", "Reference","P-value"),
       col=c("#1f78b4","grey60","#a50000"), bg="transparent", 
       lty=c(1,1,1), cex=0.8,box.lty=0,seg.len = 1,x.intersp=0.5, y.intersp=1.2)

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

par(mar=c(5.5, 5.5,4,4))
plot_quantile(as_spectra(CV_boot_Jena_PL), total_prob = 1, border = FALSE, 
              main="Bootstrapped CV of PLs vs Ref (UT-WT,EVs) in Glasshouse",ylab = "Coefficient of variance",xlab="Wavelength",
              col = "#1f78b4",cex.lab=1.4, cex.axis=1.4, cex.main=1.4, xaxs = "i")
plot(mean(as_spectra(CV_boot_Jena_PL)), col="#103c5a",lty = 2,add = TRUE)
plot(as_spectra(CV_Jena_Ref), col="grey60",lwd=2, add = TRUE)
abline(v=c(1000,1800), lty = 3)
plot_regions(as_spectra(CV_Jena_Ref),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add = TRUE)
legend(1950,0.21,legend=c("CV PLs", "Mean CV PLs","CV reference"),
       col=c("#1f78b4","#103c5a","grey60"), bg="transparent", 
       lty=c(1,2,1), cex=0.8,box.lty=0,seg.len = 1,x.intersp=0.5, y.intersp=1.2)

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

CHAL = t(as.data.frame(p_CHAL))
#CHAL = t(as.data.frame(p.adjust(p_CHAL, method = "BY"))) # all 1
colnames(CHAL) = 400:2500

HQT = t(as.data.frame(p_HQT))
#HQT = t(as.data.frame(p.adjust(p_HQT, method = "BY"))) # all 1
colnames(HQT) = 400:2500

LOX = t(as.data.frame(p_LOX))
#LOX = t(as.data.frame(p.adjust(p_LOX, method = "BY"))) # all 1
colnames(LOX) = 400:2500

RCA = t(as.data.frame(p_RCA))
a_RCA = t(as.data.frame(p.adjust(p_RCA, method = "BY"))) # min 0.4345
colnames(RCA) = 400:2500
colnames(a_RCA) = 400:2500

genogroup_o = t(as.data.frame(p_genegroup))
colnames(genogroup_o) = 400:2500
genogroup_a = t(as.data.frame(p.adjust(p_genegroup, method = "BY")))
colnames(genogroup_a) = 400:2500

plot(as_spectra(AOC), main = "TLs vs EV in Glasshouse", col = "#1f78b4", lwd=2,lty=2,
     ylab = "p-value", xlab="Wavelength",xaxs = "i",cex.lab=1.4, cex.axis=1.4, cex.main=1.5)
plot(as_spectra(LOX), col = "#a6cee3",lwd=2, lty=2,add=TRUE)
plot(as_spectra(RCA), col = "#33a02c", lwd=2,lty=2,add=TRUE)
plot(as_spectra(CHAL), col = "#fb9a99", lwd=2,lty=2,add=TRUE)
plot(as_spectra(HQT), col = "#e31a1c", lwd=2,lty=2,add=TRUE)
plot(as_spectra(CAD), col = "#fdbf6f", lwd=2,lty=2,add=TRUE)
plot(as_spectra(ACO), col = "#ff7f00", lwd=2,lty=2,add=TRUE)
plot(as_spectra(a_AOC), col = "#1f78b4",lwd=3, add=TRUE)
plot(as_spectra(a_RCA), col = "#33a02c", lwd=3,add=TRUE)
plot_regions(as_spectra(a_AOC),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add=TRUE)
abline(v=c(1000,1800), lty = 3)
abline(h=0.05, lty = 3, col="#a50000")
legend(1980,1.05,legend=c("irAOC","irAOC (adjusted)","irLOX3","irRCA","irMAX2","irCHAL","irHQT","irCAD","irACO"),
       col=c("#1f78b4","#1f78b4","#a6cee3","#33a02c","#b2df8a","#fb9a99","#e31a1c","#fdbf6f","#ff7f00"), 
       bg="transparent", lty=c(2,1,2,2,2,2,2,2,2), cex=0.8,box.lty=0,seg.len = 0.8,x.intersp=0.5,y.intersp = 1.2)

#Bootstrapping
n = 100
perm = sapply(1:n, function(i){
  boot_id = sample(nrow(Jena_TL),nrow(Jena_EV),replace=FALSE)
  boot_Jena_TL = Jena_TL[boot_id,51:2151] 
  t(as.data.frame(apply(boot_Jena_TL,2,sd) / apply(boot_Jena_TL,2,mean)))
})

CV_boot_Jena_TL = t(perm)
colnames(CV_boot_Jena_TL) = 400:2500

par(mar=c(5.5, 5.5,4, 4))
plot_quantile(as_spectra(CV_boot_Jena_TL), total_prob = 1, border = FALSE, 
              main="Bootstrapped CV of TLs vs reference (EVs) in Glasshouse",ylab = "Coefficient of variance",xlab="Wavelength",
              col = "#33a02c",cex.lab=1.4, cex.axis=1.4, cex.main=1.5, xaxs = "i")
plot(mean(as_spectra(CV_boot_Jena_TL)), col="#1a5216",lty = 2,add = TRUE)
plot(as_spectra(CV_Jena_EV), col="grey60",lwd=2, add = TRUE)
abline(v=c(1000,1800), lty = 3)
plot_regions(as_spectra(CV_Jena_EV),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add = TRUE)
legend(1950,0.03,legend=c("CV TLs", "Mean CV TLs","CV EV"),
       col=c("#33a02c","#1a5216","grey60"), bg="transparent", 
       lty=c(1,2,1), cex=0.8,box.lty=0,seg.len = 1,x.intersp=0.5, y.intersp=1.2)

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

variable = UT2[,1067-350+1]
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

noon = vector()
pm = vector()
p_ACO = vector()
p_AOC = vector()
p_CAD = vector()
p_CHAL = vector()
p_HQT = vector()
p_LOX = vector()
p_RCA = vector()

for (i in 1:2101){
  variable = UT2[,i+50]
  M = lme(variable ~ 1 + Ctime * Genotype_ID, data=UT2, random=~1|Plant_ID, method = "REML",control=lmeControl(returnObject=TRUE))
  noon[i] = summary(M)$tTable[,5][2]
  pm[i] = summary(M)$tTable[,5][3]
  p_ACO[i] = summary(M)$tTable[,5][4]
  p_AOC[i] = summary(M)$tTable[,5][5]
  p_CAD[i] = summary(M)$tTable[,5][6]
  p_CHAL[i] = summary(M)$tTable[,5][7]
  p_HQT[i] = summary(M)$tTable[,5][8]
  p_LOX[i] = summary(M)$tTable[,5][9]
  p_RCA[i] = summary(M)$tTable[,5][10]
}

# Genotypic effect
ACO = t(as.data.frame(p_ACO))
#ACO = t(as.data.frame(p.adjust(p_ACO, method = "BY"))) #all 1
colnames(ACO) = 400:2500

AOC = t(as.data.frame(p_AOC))
#AOC = t(as.data.frame(p.adjust(p_AOC, method = "BY"))) # all 1
colnames(AOC) = 400:2500

CAD = t(as.data.frame(p_CAD))
#CAD = t(as.data.frame(p.adjust(p_CAD, method = "BY"))) #all 1
colnames(CAD) = 400:2500

CHAL = t(as.data.frame(p_CHAL))
CHAL_a = t(as.data.frame(p.adjust(p_CHAL, method = "BY"))) # min 0.109, but likely due to noise
colnames(CHAL) = 400:2500
colnames(CHAL_a) = 400:2500

HQT = t(as.data.frame(p_HQT))
#HQT = t(as.data.frame(p.adjust(p_HQT, method = "BY"))) # all 1
colnames(HQT) = 400:2500

LOX = t(as.data.frame(p_LOX))
LOX_a = t(as.data.frame(p.adjust(p_LOX, method = "BY"))) # 0.1839
colnames(LOX) = 400:2500
colnames(LOX_a) = 400:2500

RCA = t(as.data.frame(p_RCA))
#RCA = t(as.data.frame(p.adjust(p_RCA, method = "BY"))) # all 1
colnames(RCA) = 400:2500

plot(as_spectra(AOC), main = "TLs vs EV in Field_UT", ylim = c(0,1.01), xaxs = "i",
     col = "#1f78b4", lty=2,lwd=2,ylab = "P-value", xlab="Wavelength",cex.lab=1.4, cex.axis=1.4, cex.main=1.5)
plot(as_spectra(LOX), col = "#a6cee3",lty=2,lwd=2, add=TRUE)
plot(as_spectra(RCA), col = "#33a02c", lty=2,lwd=2,add=TRUE)
plot(as_spectra(CHAL), col = "#fb9a99", lty=2,lwd=2,add=TRUE)
plot(as_spectra(HQT), col = "#e31a1c", lty=2,lwd=2,add=TRUE)
plot(as_spectra(CAD), col = "#fdbf6f", lty=2,lwd=2,add=TRUE)
plot(as_spectra(ACO), col = "#ff7f00", lty=2,lwd=2,add=TRUE)
plot(as_spectra(LOX_a), col = "#a6cee3",lwd=3, add=TRUE)
#plot(as_spectra(CHAL_a), col = "#fb9a99",lwd=3, add=TRUE)
plot_regions(as_spectra(AOC),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add=TRUE)
abline(v=c(1000,1800), lty = 3)
abline(h=0.05, lty = 3, col="red")
legend(1950,0.93,legend=c("irAOC","irLOX3","irLOX3 (adjusted)","irCHAL","irHQT","irCAD","irACO"),
       col=c("#1f78b4","#a6cee3","#a6cee3","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00"), 
       bg="transparent", lty=c(2,2,1,2,2,2,2,2), cex=0.8,box.lty=0,seg.len = 1,x.intersp=0.5,y.intersp = 1.2)

#Bootstrapping
n = 100
perm = sapply(1:n, function(i){
  boot_id = sample(nrow(UT_TL),nrow(UT_EV),replace=FALSE)
  boot_UT_TL = UT_TL[boot_id,51:2151] 
  t(as.data.frame(apply(boot_UT_TL,2,sd) / apply(boot_UT_TL,2,mean)))
})

CV_boot_UT_TL = t(perm)
colnames(CV_boot_UT_TL) = 400:2500

par(mar=c(5.5, 5.5,4, 4))
plot_quantile(as_spectra(CV_boot_UT_TL), total_prob = 1, border = FALSE, 
              main="Bootstrapped CV of TLs vs reference (EV) in Field_UT",ylab = "Coefficient of variance",xlab="Wavelength",
              col = "#33a02c",cex.lab=1.4, cex.axis=1.4, cex.main=1.5, xaxs = "i")
plot(mean(as_spectra(CV_boot_UT_TL)), col="#1a5216",lty = 2,add = TRUE)
plot(as_spectra(CV_UT_EV), col="grey60",lwd=2, add = TRUE)
abline(v=c(1000,1800), lty = 3)
plot_regions(as_spectra(CV_UT_EV),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add = TRUE)
legend(1950,0.42,legend=c("CV TLs", "Mean CV TLs","CV EV"),
       col=c("#33a02c","#1a5216","grey60"), bg="transparent", 
       lty=c(1,2,1), cex=0.8,box.lty=0,seg.len = 1,x.intersp=0.5, y.intersp=1.2)

# Time effect
o_noon = t(as.data.frame(noon))
p_noon = t(as.data.frame(p.adjust(noon, method = "BY")))
colnames(o_noon) = 400:2500
colnames(p_noon) = 400:2500

o_pm = t(as.data.frame(pm))
p_pm = t(as.data.frame(p.adjust(pm, method = "BY")))
colnames(o_pm) = 400:2500
colnames(p_pm) = 400:2500
rownames(t(p_pm))[t(p_pm)< 1] #671~700,1001~2500
rownames(t(o_pm))[t(o_pm)< 0.05] #1153-1165,1266-2168

plot_quantile(as_spectra(filter(UT2,Ctime=="am")[,51:2151]),ylim = c(0,1.01), total_prob = 1, border=FALSE,col = rgb(0.992, 0.749, 0.435,0.5),main = "Effect of time in Field_UT", 
              ylab = "Reflectance", xlab="Wavelength",xaxs = "i",cex.lab=1.4, cex.axis=1.4, cex.main=1.5)
plot_quantile(as_spectra(filter(UT2,Ctime=="noon")[,51:2151]),col = rgb(1, 0.498, 0,0.5), total_prob = 1, border=FALSE,add=TRUE)
plot_quantile(as_spectra(filter(UT2,Ctime=="pm")[,51:2151]),col = rgb(0.122, 0.471, 0.706,0.5), total_prob = 1, border=FALSE,add=TRUE)
plot(as_spectra(o_noon), col = "#ff7f00", lwd=2, lty=2,add=TRUE)
plot(as_spectra(p_noon), col = "#ff7f00", lwd=2,add=TRUE)
plot(as_spectra(o_pm),col = "#1f78b4", lty=2, lwd=2,add=TRUE)
plot(as_spectra(p_pm), col = "#1f78b4",lwd=2, add=TRUE)
plot_regions(as_spectra(o_noon),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add=TRUE)
abline(v=c(1000,1800), lty = 3)
abline(h=0.05, lty = 3, col="red")
axis(side = 4, cex.axis=1.5)
mtext("p-value",side=4,line=3,cex = 1.5) 
legend(2000,1,legend=c("AM","Noon vs AM","PM vs AM"),
       col=c("#fdbf6f","#ff7f00","#1f78b4"), 
       bg="transparent", lty=c(1,1,1), cex=0.8,box.lty=0,seg.len = 0.8,x.intersp=0.5,y.intersp = 1.2)

# Model 2: UT_dev
# UT_dev leaf position
utdev = mutate(UT_dev,leaf = ifelse(Leaf_ID == 2,"a2",paste0("b",Leaf_ID)),
               Plant_ID = paste0("p",Plant_ID))

# Decide main effects on Wavelength 1399,689,717,711
variable = utdev[,1399-350+1]
M1 <- lme(variable ~ 1 + Ctime + leaf, data=utdev, random=~1|Plant_ID, method = "REML")
M2 <- lme(variable ~ 1 + Ctime * leaf, data=utdev, random=~1|Plant_ID, method = "REML")
AIC(M1,M2) #M1 < M2

variable = utdev[,689-350+1]
M1 <- lme(variable ~ 1 + Ctime + leaf, data=utdev, random=~1|Plant_ID, method = "REML")
M2 <- lme(variable ~ 1 + Ctime * leaf, data=utdev, random=~1|Plant_ID, method = "REML")
AIC(M1,M2) #M1 < M2 

variable = utdev[,717-350+1]
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
p_leaf1 = vector()
p_leaf3 = vector()
avo_leaf = vector()

for (i in 1:2101){
  variable = utdev[,i+50]
  M = lme(variable ~ 1 + Ctime + leaf, data=utdev, random=~1|Plant_ID, method = "REML")
  p_noon[i] = summary(M)$tTable[,5][2]
  p_leaf1[i] = summary(M)$tTable[,5][3]
  p_leaf3[i] = summary(M)$tTable[,5][4]
  
  A = anova(M)
  avo_leaf[i] = A$`p-value`[3]
}

leaf1_o = t(as.data.frame(p_leaf1))
leaf1_a = t(as.data.frame(p.adjust(p_leaf1, method = "BY")))
colnames(leaf1_o) = 400:2500
colnames(leaf1_a) = 400:2500

leaf3_o = t(as.data.frame(p_leaf3))
leaf3_a = t(as.data.frame(p.adjust(p_leaf3, method = "BY")))
colnames(leaf3_o) = 400:2500
colnames(leaf3_a) = 400:2500
sigwavs(leaf3_o) #545-711


avo_leaf_o = t(as.data.frame(avo_leaf))
avo_leaf_a = t(as.data.frame(p.adjust(avo_leaf, method = "BY")))
colnames(avo_leaf_o) = 400:2500
sigwavs(avo_leaf_o) #488-722
colnames(avo_leaf_a) = 400:2500 #min 0.2

plot_quantile(as_spectra(filter(utdev,leaf=="b1")[,51:2151]),ylim = c(0,1.01), total_prob = 1, border=FALSE,col = rgb(0.992, 0.749, 0.435,0.5),main = "Effects of leaf position in Field_UT", 
              xaxs = "i",ylab = "Reflectance", xlab="Wavelength",cex.lab=1.4, cex.axis=1.4, cex.main=1.5)
plot_quantile(as_spectra(filter(utdev,leaf=="a2")[,51:2151]),col = rgb(1, 0.498, 0,0.5), total_prob = 1, border=FALSE,add=TRUE)
plot_quantile(as_spectra(filter(utdev,leaf=="b3")[,51:2151]),col = rgb(0.122, 0.471, 0.706,0.5), total_prob = 1, border=FALSE,add=TRUE)
plot(as_spectra(leaf1_o), col = "#ff7f00", lwd = 2,lty=2, add=TRUE)
plot(as_spectra(leaf1_a), col = "#ff7f00", lwd=2, add=TRUE)
plot(as_spectra(leaf3_o),col = "#1f78b4", lwd=2, lty=2, add=TRUE)
plot(as_spectra(leaf3_a), col = "#1f78b4", lwd=2, add=TRUE)
abline(v=c(1000,1800), lty = 3)
abline(h=0.05, lty = 3, col="red")
plot_regions(as_spectra(leaf1_o),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add=TRUE)
axis(side = 4, cex.axis=1.5)
mtext("p-value",side=4,line=3,cex = 1.5) 

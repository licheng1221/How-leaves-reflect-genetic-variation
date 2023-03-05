library(tidyverse)
library(spectrolab)
library(factoextra)

##### Arizona ##### 
A = readRDS("~/Desktop/Nicotiana/Final/Data_rds/A_AZ.rds")
R = readRDS("~/Desktop/Nicotiana/Final/Data_rds/R_AZ.rds")

## Coefficient of variation
SV_all = 100* apply(R[51:2151],2,sd) / apply(R[51:2151],2,mean)
SV_all = t(as.data.frame(SV_all))
colnames(SV_all) = 400:2500

SV_ref = 100* apply(filter(R,Genotype_ID == "UT-WT")[51:2151],2,sd) / apply(filter(R,Genotype_ID == "UT-WT")[51:2151],2,mean)
SV_ref = t(as.data.frame(SV_ref))
colnames(SV_ref) = 400:2500

SV_ril = 100* apply(filter(R,startsWith(Genotype_ID,"M"))[51:2151],2,sd) / apply(filter(R,startsWith(Genotype_ID,"M"))[51:2151],2,mean)
SV_ril = t(as.data.frame(SV_ril))
colnames(SV_ril) = 400:2500

SV_pl = 100* apply(filter(R,startsWith(Genotype_ID,"P"))[51:2151],2,sd) / apply(filter(R,startsWith(Genotype_ID,"P"))[51:2151],2,mean)
SV_pl = t(as.data.frame(SV_pl))
colnames(SV_pl) = 400:2500

# Absolute uncertainty
mean_AU = t(apply(A[,51:2151],2,mean))
colnames(mean_AU) = 400:2500

# Relative uncertainty
RU = 100 * A[,51:2151] / R[51:2151]
mean_RU = apply(RU, 2, mean)
mean_RU = t(as.data.frame(mean_RU))
colnames(mean_RU) = 400:2500

# Plots
plot(as_spectra(SV_all), main="Coefficient of variantion Field_AZ",ylab = "Cofficient of Variation (%)", xlab = "Wavelength",
     ylim = c(0,20),col = "#00AFBB",xaxs = "i", yaxs = "i",lwd=2,cex.lab=1.4, cex.axis=1.4, cex.main=1.5)
plot(as_spectra(SV_ref),col = "grey60", lwd=2, add=TRUE)
plot(as_spectra(SV_ril),col = "#a6cee3", lwd=2,add=TRUE)
plot(as_spectra(SV_pl),col = "#1f78b4", lwd=2,add=TRUE)
plot(as_spectra(mean_RU),col = "#9400D3",lwd=2,add=TRUE)
abline(v=c(1000,1800), lty = 3)
plot_regions(as_spectra(SV_all),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2))
legend(1750,5,legend=c("SV All","SV PLs","SV RILs","SV Ref","Absolute Uncertainty"),
       col=c("#00AFBB","#1f78b4","#a6cee3","grey60","#9400D3"), bg="transparent", 
       lty=1, cex=0.8,box.lty=0,seg.len = 1,x.intersp=0.5, y.intersp=1.2)

## PCA
AZ = mutate(R, Genos = ifelse(Genotype_ID == "UT-WT", "Ref", "non-Ref"),
            GeneGroup = ifelse(Genotype_ID == "UT-WT", "Ref", 
                               ifelse(startsWith(Genotype_ID,"P"),"PL","RIL")))

pca.AZ = prcomp(AZ[,51:2151],scale. = TRUE)
fviz_pca_ind(pca.AZ, axes = c(1,2),label = "none", habillage = AZ$Genos, palette = c("#00AFBB", "grey60"),addEllipses = TRUE) + 
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(.15, .85),
        text=element_text(size=20),
        legend.text = element_text(size=15),
        legend.background = element_rect(fill="transparent"))
fviz_pca_ind(pca.AZ, axes = c(3,4),label = "none", habillage = AZ$Genos, palette = c("#00AFBB", "grey60"),addEllipses = TRUE) + 
        theme_classic()+
        theme(plot.title = element_text(hjust = 0.5),
              legend.position = c(.15, .85),
              text=element_text(size=20),
              legend.text = element_text(size=15),
              legend.background = element_rect(fill="transparent"))

fviz_eig(pca.AZ, addlabels = TRUE, main="Field_AZ", barfill ="#00AFBB", barcolor = "#00AFBB", ncp = 8) + 
        ylim(0,70) + 
        theme(plot.title = element_text(hjust = 0.5),
              text = element_text(size = 20))
which.max((get_pca_var(pca.AZ)$contrib)[,1]) # 1877, 1478
which.max((get_pca_var(pca.AZ)$contrib)[,2]) # 824, 425
which.max((get_pca_var(pca.AZ)$contrib)[,3]) # 693, 294
which.max((get_pca_var(pca.AZ)$contrib)[,4]) # 552, 153

##### Jena #####
A_2 = readRDS("~/Desktop/Nicotiana/Final/Data_rds/A_Jena.rds")
R_2 = readRDS("~/Desktop/Nicotiana/Final/Data_rds/R_Jena.rds")

## Coefficient of variation
SV_all_2 = 100* apply(R_2[51:2151],2,sd) / apply(R_2[51:2151],2,mean)
SV_all_2 = t(as.data.frame(SV_all_2))
colnames(SV_all_2) = 400:2500

SV_ref_2 = 100* apply(filter(R_2,Genotype_ID %in% c("UT_WT","pRESC2NC","pSOL3NC"))[51:2151],2,sd) / apply(filter(R_2,Genotype_ID %in% c("UT_WT","pRESC2NC","pSOL3NC"))[51:2151],2,mean)
SV_ref_2 = t(as.data.frame(SV_ref_2))
colnames(SV_ref_2) = 400:2500

SV_TL_2 = 100* apply(filter(R_2,startsWith(Genotype_ID,"ir"))[51:2151],2,sd) / apply(filter(R_2,startsWith(Genotype_ID,"ir"))[51:2151],2,mean)
SV_TL_2 = t(as.data.frame(SV_TL_2))
colnames(SV_TL_2) = 400:2500

SV_PL_2 = 100* apply(filter(R_2,startsWith(Genotype_ID,"P"))[51:2151],2,sd) / apply(filter(R_2,startsWith(Genotype_ID,"P"))[51:2151],2,mean)
SV_PL_2 = t(as.data.frame(SV_PL_2))
colnames(SV_PL_2) = 400:2500

# Absolute uncertainty
mean_AU_2 = t(apply(A_2[,51:2151],2,mean))
colnames(mean_AU_2) = 400:2500

# Relative uncertainty
RU_2 = 100 * A_2[,51:2151] / R_2[51:2151]
mean_RU_2 = apply(RU_2, 2, mean)
mean_RU_2 = t(as.data.frame(mean_RU_2))
colnames(mean_RU_2) = 400:2500

plot(as_spectra(SV_all_2), main="Coefficient of variation of samples in Glasshouse",col = "#E7B800",
     ylab = "Coefficient of Variation (%)", xlab = "Wavelength",ylim = c(0,20),lwd=2, xaxs = "i", yaxs = "i",cex.lab=1.4, cex.axis=1.4, cex.main=1.5)
plot(as_spectra(SV_ref_2),col = "grey60", lwd=2,add=TRUE)
plot(as_spectra(SV_TL_2),col = "#33a02c",lwd=2, add=TRUE)
plot(as_spectra(SV_PL_2),col = "#1f78b4", lwd=2,add=TRUE)
plot(as_spectra(mean_RU_2), col = "#9400D3",lwd=2,add=TRUE)
abline(v=c(1000,1800), lty = 3)
plot_regions(as_spectra(SV_all_2),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add = TRUE)
legend(1750,5,legend=c("All","TLs","PLs","Ref","Absolute Uncertainty"),
       col=c("#E7B800","#33a02c","#1f78b4","grey60","#9400D3"), bg="transparent", 
       lty=1, cex=0.8,box.lty=0,seg.len = 1,x.intersp=0.5,y.intersp = 1.2)

## PCA
Jena = mutate(R_2, Genos = ifelse(Genotype_ID %in% c("UT_WT","pRESC2NC","pSOL3NC"), "Ref", "non-Ref"),
              GeneGroup = ifelse(Genos == "Ref", "Ref", 
                                 ifelse(startsWith(Genotype_ID,"ir"),"TL","PL")))

pca.Jena = prcomp(Jena[,51:2151],scale. = TRUE)
fviz_pca_ind(pca.Jena, axes = c(1,2),label = "none", habillage = Jena$Genos, palette = c("#E7B800", "grey60"),addEllipses = TRUE) + 
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(.15, .85),
        text=element_text(size=20),
        legend.text = element_text(size=15),
        legend.background = element_rect(fill="transparent"))

fviz_pca_ind(pca.Jena, axes = c(3,4),label = "none", habillage = Jena$Genos, palette = c("#E7B800", "grey60"),addEllipses = TRUE) + 
        theme_classic()+
        theme(plot.title = element_text(hjust = 0.5),
              legend.position = c(.15, .85),
              text=element_text(size=20),
              legend.text = element_text(size=15),
              legend.background = element_rect(fill="transparent"))

fviz_eig(pca.Jena, addlabels = TRUE, main="Glasshouse", barfill ="#E7B800", barcolor = "#E7B800", ncp = 8) + 
        ylim(0,70) + 
        theme(plot.title = element_text(hjust = 0.5),
              text = element_text(size = 20))

which.max((get_pca_var(pca.Jena)$contrib)[,1]) # 1354, 955
which.max((get_pca_var(pca.Jena)$contrib)[,2]) # 731, 332
which.max((get_pca_var(pca.Jena)$contrib)[,3]) # 654, 255
which.max((get_pca_var(pca.Jena)$contrib)[,4]) # 441, 42

##### Utah #####
A_3 = readRDS("~/Desktop/Nicotiana/Final/Data_rds/A_UT_lin.rds")
R_3 = readRDS("~/Desktop/Nicotiana/Final/Data_rds/R_UT_lin.rds")

## Coefficient of variation
SV_all_3 = 100* apply(R_3[51:2151],2,sd) / apply(R_3[51:2151],2,mean)
SV_all_3 = t(as.data.frame(SV_all_3))
colnames(SV_all_3) = 400:2500

SV_ref_3 = 100* apply(filter(R_3,Genotype_ID %in% c("pRESC2NC"))[51:2151],2,sd) / apply(filter(R_3,Genotype_ID %in% c("pRESC2NC"))[51:2151],2,mean)
SV_ref_3 = t(as.data.frame(SV_ref_3))
colnames(SV_ref_3) = 400:2500

SV_TL_3 = 100* apply(filter(R_3,startsWith(Genotype_ID,"ir"))[51:2151],2,sd) / apply(filter(R_3,startsWith(Genotype_ID,"ir"))[51:2151],2,mean)
SV_TL_3 = t(as.data.frame(SV_TL_3))
colnames(SV_TL_3) = 400:2500

SV_all_am = 100* apply(filter(R_3,Ctime=="am")[,51:2151],2,sd) / apply(filter(R_3,Ctime=="am")[,51:2151],2,mean)
SV_all_am = t(as.data.frame(SV_all_am))
colnames(SV_all_am) = 400:2500

SV_all_noon = 100* apply(filter(R_3,Ctime=="noon")[,51:2151],2,sd) / apply(filter(R_3,Ctime=="noon")[,51:2151],2,mean)
SV_all_noon = t(as.data.frame(SV_all_noon))
colnames(SV_all_noon) = 400:2500

SV_all_pm = 100* apply(filter(R_3,Ctime=="pm")[,51:2151],2,sd) / apply(filter(R_3,Ctime=="pm")[,51:2151],2,mean)
SV_all_pm = t(as.data.frame(SV_all_pm))
colnames(SV_all_pm) = 400:2500

# Absolute uncertainty
mean_AU_3 = t(apply(A_3[,51:2151],2,mean))
colnames(mean_AU_3) = 400:2500

# Relative uncertainty
RU_3 = 100 * A_3[,51:2151] / R_3[51:2151]
mean_RU_3 = apply(RU_3, 2, mean)
mean_RU_3 = t(as.data.frame(mean_RU_3))
colnames(mean_RU_3) = 400:2500

plot(as_spectra(SV_all_3), main="Coefficient of variation of samples in Field_UT",ylab = "Coefficient of Variation (%)", xlab = "Wavelength",
     ylim = c(0,40), col = "#FC4E07", lwd=2,xaxs = "i", yaxs = "i",cex.lab=1.4, cex.axis=1.4, cex.main=1.5)
#areaplot(t(SV_all_3),col = rgb(252,78,7,max=255,alpha = 100),add=TRUE)
plot(as_spectra(SV_ref_3),col = "grey60", lwd=2,add=TRUE)
plot(as_spectra(SV_TL_3),col = "#33a02c", lwd=2,add=TRUE)
plot(as_spectra(SV_all_am), col = "#fdbf6f",lwd=2,add=TRUE)
plot(as_spectra(SV_all_noon), col = "#ff7f00",lwd=2,add=TRUE)
plot(as_spectra(SV_all_pm), col = "#1f78b4",lwd=2,add=TRUE)
plot(as_spectra(mean_RU_3), col = "#9400D3",lwd=2,add=TRUE)
abline(v=c(1000,1800), lty = 3)
plot_regions(as_spectra(SV_all_2),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add = TRUE)
legend(1750,80,legend=c("All","TLs","Ref","am","noon","pm","Absolute Uncertainty"),
       col=c("#FC4E07","#33a02c","grey60","#fdbf6f","#ff7f00","#1f78b4","#9400D3"), bg="transparent", 
       lty=1, cex=0.8,box.lty=0,seg.len = 1,x.intersp=0.5,y.intersp = 1.2)

## PCA
R_4 = readRDS("~/Desktop/Nicotiana/Final/Data_rds/R_UT_dev.rds")
# remove colume Leaf_ID for now
R_4 = R_4[,-2154]
UT = rbind(R_3,R_4) %>%
        mutate(Genos = ifelse(Genotype_ID %in% c("UT_WT","pRESC2NC","pSOL3NC"), "Ref", "non-Ref"),
               GeneGroup = ifelse(Genos == "Ref","Ref","TL"))

pca.UT = prcomp(UT[,51:2151],scale. = TRUE)
fviz_pca_ind(pca.UT, axes = c(1,2),label = "none", habillage = UT$Genos, palette = c("#FC4E07", "grey60"),addEllipses = TRUE) + 
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(.15, .85),
        text=element_text(size=20),
        legend.text = element_text(size=15),
        legend.background = element_rect(fill="transparent"))
fviz_pca_ind(pca.UT, axes = c(3,4),label = "none", habillage = UT$Genos, palette = c("#FC4E07", "grey60"),addEllipses = TRUE) + 
        theme_classic()+
        theme(plot.title = element_text(hjust = 0.5),
              legend.position = c(.15, .85),
              text=element_text(size=20),
              legend.text = element_text(size=15),
              legend.background = element_rect(fill="transparent"))

fviz_eig(pca.UT, addlabels = TRUE, main="Field_UT", barfill ="#FC4E07", barcolor = "#FC4E07", ncp = 8) + 
        ylim(0,75) + 
        theme(plot.title = element_text(hjust = 0.5),
              text = element_text(size =20))



# PCA UT_lin
pca.UT.lin = prcomp(R_3[,51:2151],scale. = TRUE)

which.max((get_pca_var(pca.UT.lin)$contrib)[,1]) # 1516,1117
which.max((get_pca_var(pca.UT.lin)$contrib)[,2]) # 1005,606
which.max((get_pca_var(pca.UT.lin)$contrib)[,3]) # 536,137
which.max((get_pca_var(pca.UT.lin)$contrib)[,4]) # 1801,1402

# PCA UT_dev
pca.UT.dev = prcomp(R_4[,51:2151],scale. = TRUE)

which.max((get_pca_var(pca.UT.dev)$contrib)[,1]) # 1399,1000
which.max((get_pca_var(pca.UT.dev)$contrib)[,2]) # 689,290
which.max((get_pca_var(pca.UT.dev)$contrib)[,3]) # 717,318
which.max((get_pca_var(pca.UT.dev)$contrib)[,4]) # 711,312

# Measurement of uncertainty comparision
#relative
plot(as_spectra(mean_RU_3), main="Measurement uncertainty in three envs",col = "#FC4E07",
     ylab = "Measurement of uncertainty (%)", xlab = "Wavelength",ylim = c(0,1.1),lwd=2, xaxs = "i", yaxs = "i",cex.lab=1.4, cex.axis=1.4, cex.main=1.5)
plot(as_spectra(mean_RU),col = "#00AFBB", lwd=2,add=TRUE)
plot(as_spectra(mean_RU_2),col = "#E7B800",lwd=2, add=TRUE)
abline(v=c(1000,1800), lty = 3)
plot_regions(as_spectra(mean_RU),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add = TRUE)
#absolute
plot(as_spectra(mean_AU_3), main="Measurement uncertainty in three envs",col = "#FC4E07",
     ylab = "Measurement of uncertainty (%)", xlab = "Wavelength",ylim = c(0,0.003),lwd=2, xaxs = "i", yaxs = "i",cex.lab=1.4, cex.axis=1.4, cex.main=1.5)
plot(as_spectra(mean_AU),col = "#00AFBB", lwd=2,add=TRUE)
plot(as_spectra(mean_AU_2),col = "#E7B800",lwd=2, add=TRUE)
abline(v=c(1000,1800), lty = 3)
plot_regions(as_spectra(mean_RU),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add = TRUE)

## Contributions to PC 1~3
spec_AZ = as_spectra(t(get_pca_var(pca.AZ)$contrib))
spec_Jena = as_spectra(t(get_pca_var(pca.Jena)$contrib))
spec_UT = as_spectra(t(get_pca_var(pca.UT)$contrib))

plot(spec_Jena[1,], 
     main="Contribution of variables to Dimension 1",
     xlab="Wavelength", ylab = "Contribution(%)",xaxs = "i", 
     col = "#E7B800",lwd=2,cex.lab=1.4, cex.axis=1.4, cex.main=1.5)
plot(spec_AZ[1,],lwd=2,col="#00AFBB",add=TRUE)
plot(spec_UT[1,],lwd=2,col="#FC4E07",add=TRUE)
plot_regions(spec_Jena, col = grDevices::rgb(0.7, 0.7, 0.7, 0.2),add=TRUE)
abline(v=c(1000,1800), lty = 3)
legend(1950,0.01,legend=c("Field_AZ", "Field_UT","Glasshouse"),
       col=c("#00AFBB","#FC4E07","#E7B800"), bg="transparent", lty=1, cex=0.8,box.lty=0,seg.len = 1,x.intersp=0.5, y.intersp=1)

plot(spec_AZ[2,], 
     main="Contribution of variables to Dimension 2",
     xlab="Wavelength", ylab = "Contribution(%)",xaxs = "i", 
     col = "#00AFBB",lwd=2,cex.lab=1.4, cex.axis=1.4, cex.main=1.5)
plot(spec_Jena[2,],lwd=2,col="#E7B800",add=TRUE)
plot(spec_UT[2,],lwd=2,col="#FC4E07",add=TRUE)
plot_regions(spec_Jena, col = grDevices::rgb(0.7, 0.7, 0.7, 0.2),add=TRUE)
abline(v=c(1000,1800), lty = 3)
legend(1900,0.17,legend=c("Field_AZ", "Field_UT","Glasshouse"),
       col=c("#00AFBB","#FC4E07","#E7B800"), bg="transparent", lty=1, cex=0.8,box.lty=0,seg.len = 1,x.intersp=0.5, y.intersp=1)

plot(spec_Jena[3,],col = "#E7B800",
     main="Contribution of variables to Dimension 3",
     xlab="Wavelength", ylab = "Contribution(%)",xaxs = "i", 
     cex.lab=1.4, cex.axis=1.4,lwd=2, cex.main=1.5)
plot(spec_AZ[3,],col="#00AFBB",lwd=2,add=TRUE)
plot(spec_UT[3,],col="#FC4E07",lwd=2,add=TRUE)
plot_regions(spec_Jena, col = grDevices::rgb(0.7, 0.7, 0.7, 0.2),add=TRUE)
abline(v=c(1000,1800), lty = 3)
legend(1900,0.34,legend=c("Field_AZ", "Field_UT","Glasshouse"),
       col=c("#00AFBB","#FC4E07","#E7B800"), bg="transparent", lty=1, cex=0.8,box.lty=0,seg.len = 1,x.intersp=0.5, y.intersp=1)

plot(spec_AZ[4,],col="#00AFBB",
     main="Contribution of variables to Dimension 4",
     xlab="Wavelength", ylab = "Contribution(%)",xaxs = "i", 
     cex.lab=1.4, cex.axis=1.4,lwd=2, cex.main=1.5)
plot(spec_UT[4,],col="#FC4E07",lwd=2,add=TRUE)
plot(spec_Jena[4,],col = "#E7B800",lwd=2,add=TRUE)
plot_regions(spec_Jena, col = grDevices::rgb(0.7, 0.7, 0.7, 0.2),add=TRUE)
abline(v=c(1000,1800), lty = 3)
legend(1900,0.34,legend=c("Field_AZ", "Field_UT","Glasshouse"),
       col=c("#00AFBB","#FC4E07","#E7B800"), bg="transparent", lty=1, cex=0.8,box.lty=0,seg.len = 1,x.intersp=0.5, y.intersp=1)

## All PCA
All = as.data.frame(rbind(AZ[,c(51:2152,2158:2160)],Jena[,c(51:2152,2156:2158)],UT[,c(51:2152,2155:2157)]))
pca.all = prcomp(All[,1:2101],scale. = TRUE)

fviz_pca_ind(pca.all, label = "none", title = "PCA - All samples", habillage = All$Site, palette = c("#00AFBB", "#FC4E07","#E7B800"),addEllipses = TRUE) + 
        theme_classic()+
        theme(plot.title = element_text(hjust = 0.5),
              legend.position = c(.85, .92),
              text=element_text(size=20),
              legend.text = element_text(size=15),
              legend.background = element_rect(fill="transparent"))
fviz_pca_ind(pca.all, axes = c(3,4),label = "none", title = "PCA - All samples", habillage = All$Site, palette = c("#00AFBB", "#FC4E07","#E7B800"),addEllipses = TRUE) + 
        theme_classic()+
        theme(plot.title = element_text(hjust = 0.5),
              legend.position = c(.85, .92),
              text=element_text(size=20),
              legend.text = element_text(size=15),
              legend.background = element_rect(fill="transparent"))

# PCA UT_WT
All_UT = filter(All,Genotype_ID %in% c("UT-WT","UT_WT"))
pca.all.ut = prcomp(All_UT[,1:2101],scale. = TRUE)
fviz_pca_ind(pca.all.ut, label = "none", title = "PCA - UT-WT samples",habillage = All_UT$Site, palette = c("#00AFBB", "#FC4E07","#E7B800"), addEllipses = TRUE) + 
        theme_classic()+
        theme(plot.title = element_text(hjust = 0.5),
              legend.position = c(.85, .15),
              text=element_text(size=20),
              legend.text = element_text(size=15),
              legend.background = element_rect(fill="transparent"))
fviz_pca_ind(pca.all.ut, axes = c(3,4),label = "none", title = "PCA - UT-WT samples",habillage = All_UT$Site, palette = c("#00AFBB", "#FC4E07","#E7B800"), addEllipses = TRUE) + 
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(.85, .15),
        text=element_text(size=20),
        legend.text = element_text(size=15),
        legend.background = element_rect(fill="transparent"))

##PCs, contribution of variable to dimension 1~4, PCA_All
fviz_eig(pca.all, addlabels = TRUE, main="All experiments", barfill ="grey60", barcolor = "grey60", ncp = 8) + 
        ylim(0,70) + 
        theme(plot.title = element_text(hjust = 0.5),
              text = element_text(size =20))

spec_all = as_spectra(t(get_pca_var(pca.all)$contrib))
plot(spec_all[4,], 
     main="Contribution of variables to dimension 1~4",xaxs = "i",
     xlab="Wavelength", ylab = "Contribution(%)",
     col = "grey60",lwd=2,cex.lab=1.4, cex.axis=1.4, cex.main=1.5)
plot(spec_all[1,],lwd=2,col="#00AFBB",add=TRUE)
plot(spec_all[2,],lwd=2,col="#FC4E07",add=TRUE)
plot(spec_all[3,],lwd=2,col="#E7B800",add=TRUE)
plot_regions(spec_all[1,], col = grDevices::rgb(0.7, 0.7, 0.7, 0.2),add=TRUE)
legend(2000,0.18,legend=c("Dimension 1", "Dimension 2","Dimension 3","Dimension 4"),
       col=c("#00AFBB","#FC4E07","#E7B800","grey60"), bg="transparent", lty=1, cex=0.8,box.lty=0,seg.len = 1,x.intersp=0.5, y.intersp=0.8)

# PCs, ontribution of variable to dimension 1~4, PCA_UT-WT
fviz_eig(pca.all.ut, addlabels = TRUE, main="All experiments", barfill ="grey60", barcolor = "grey60", ncp = 8) + 
  ylim(0,70) + 
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size =20))

spec_all_ut = as_spectra(t(get_pca_var(pca.all.ut)$contrib))
plot(spec_all_ut[4,], 
     main="Contribution of variables to dimension 1~4",
     xlab="Wavelength", ylab = "Contribution(%)",xaxs = "i",
     col = "grey60",lwd=2,cex.lab=1.4, cex.axis=1.4, cex.main=1.5)
plot(spec_all_ut[1,],lwd=2,col="#00AFBB",add=TRUE)
plot(spec_all_ut[2,],lwd=2,col="#FC4E07",add=TRUE)
plot(spec_all_ut[3,],lwd=2,col="#E7B800",add=TRUE)
plot_regions(spec_all_ut[1,], col = grDevices::rgb(0.7, 0.7, 0.7, 0.2),add=TRUE)
legend(2000,0.18,legend=c("Dimension 1", "Dimension 2","Dimension 3","Dimension 4"),
       col=c("#00AFBB","#FC4E07","#E7B800","grey60"), bg="transparent", lty=1, cex=0.8,box.lty=0,seg.len = 1,x.intersp=0.5, y.intersp=0.8)

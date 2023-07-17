library(tidyverse)
library(spectrolab)
library(factoextra)
library(cowplot)
library(gridGraphics)
#library(areaplot)

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

# Fig 5a
plot(as_spectra(SV_all), main="CV Field_AZ",ylab = "CV (%)", xlab = "Wavelength (nm)",
     ylim = c(0,20),col = "#56B4E9",xaxs = "i", yaxs = "i",lwd=3)
plot(as_spectra(SV_ref),col = "#000000", lwd=3, lty=2, add=TRUE)
plot(as_spectra(SV_ril),col = "#F0E442", lwd=3,lty=2,add=TRUE)
plot(as_spectra(SV_pl),col = "#E69F00", lwd=3,lty=2,add=TRUE)
plot(as_spectra(mean_RU),col = "#CC79A7",lwd=3,lty=3,add=TRUE) 
abline(v=c(1000,1800), lty = 3)
plot_regions(as_spectra(SV_all),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2))
legend(1900,7,legend=c("All","PLs","RILs","Ref","MU"),
       col=c("#56B4E9","#E69F00","#F0E442","#000000","#CC79A7"), bg="transparent", 
       lty=c(1,2,2,2,3), lwd=3,cex=1,box.lty=0,seg.len = 2,x.intersp=0.5, y.intersp=0.5)
Fig_5a <- recordPlot()

## PCA
AZ = mutate(R, Genos = ifelse(Genotype_ID == "UT-WT", "Ref", "non-Ref"),
            GeneGroup = ifelse(Genotype_ID == "UT-WT", "Ref", 
                               ifelse(startsWith(Genotype_ID,"P"),"PL","RIL")))

# Fig S2d
pca.AZ = prcomp(AZ[,51:2151],scale. = TRUE)
Fig_S2d <- fviz_pca_ind(pca.AZ, axes = c(1,2),label = "none", habillage = AZ$Genos, palette = c("#56B4E9", "#000000"),addEllipses = TRUE) + 
  my_theme +
  theme(legend.position = c(.11, .90))+
  ggtitle("PCA Field_AZ")

#Fig 5e
Fig_5e <- fviz_pca_ind(pca.AZ, axes = c(3,4),label = "none", habillage = AZ$Genos, palette = c("#56B4E9", "#000000"),addEllipses = TRUE) + 
  my_theme +
  theme(legend.position = c(.11, .90)) +
  ggtitle("PCA Field_AZ")

#Fig S2a
Fig_S2a <- fviz_eig(pca.AZ, addlabels = TRUE, main="Field_AZ", barfill ="#56B4E9", barcolor = "#56B4E9", ncp = 8) + 
  ylim(0,70) + 
  xlab("PCs") +
  my_theme +
  ggtitle("Explained variance Field_AZ")

which.max((get_pca_var(pca.AZ)$contrib)[,1]) # 1877
which.max((get_pca_var(pca.AZ)$contrib)[,2]) # 824
which.max((get_pca_var(pca.AZ)$contrib)[,3]) # 693
which.max((get_pca_var(pca.AZ)$contrib)[,4]) # 552

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

#Fig 5c
plot(as_spectra(SV_all_2), main="CV Glasshouse",col = "#009E73",
     ylab = "CV (%)", xlab = "Wavelength (nm)",ylim = c(0,20),lwd=3, xaxs = "i", yaxs = "i")
plot(as_spectra(SV_ref_2),col = "#000000", lwd=3,lty=2,add=TRUE)
plot(as_spectra(SV_TL_2),col = "#0072B2",lwd=3,lty=2, add=TRUE)
plot(as_spectra(SV_PL_2),col = "#E69F00", lwd=3,lty=2,add=TRUE)
plot(as_spectra(mean_RU_2), col = "#CC79A7",lwd=3,lty=3,add=TRUE)
abline(v=c(1000,1800), lty = 3)
plot_regions(as_spectra(SV_all_2),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add = TRUE)
legend(1900,7,legend=c("All","TLs","PLs","Ref","MU"),
       col=c("#009E73","#0072B2","#E69F00","#000000","#CC79A7"), bg="transparent", 
       lty=c(1,2,2,2,3), lwd=3,cex=1,box.lty=0,seg.len = 2,x.intersp=0.5,y.intersp = 0.5)
Fig_5c <- recordPlot()

## PCA
Jena = mutate(R_2, Genos = ifelse(Genotype_ID %in% c("UT_WT","pRESC2NC","pSOL3NC"), "Ref", "non-Ref"),
              GeneGroup = ifelse(Genos == "Ref", "Ref", 
                                 ifelse(startsWith(Genotype_ID,"ir"),"TL","PL")))

pca.Jena = prcomp(Jena[,51:2151],scale. = TRUE)
# Fig 5g
Fig_5g <- fviz_pca_ind(pca.Jena, axes = c(1,2),label = "none", habillage = Jena$Genos, palette = c("#009E73", "#000000"),addEllipses = TRUE) + 
  my_theme +
  theme(legend.position = c(.11, .90)) +
  ggtitle("PCA Glasshouse")
# Fig S2f
Fig_S2f <- fviz_pca_ind(pca.Jena, axes = c(3,4),label = "none", habillage = Jena$Genos, palette = c("#009E73", "#000000"),addEllipses = TRUE) + 
  my_theme +
  theme(legend.position = c(.11, .90))+
  ggtitle("PCA Glasshouse")
# Fig S2c
Fig_S2c <- fviz_eig(pca.Jena, addlabels = TRUE, main="Glasshouse", barfill ="#009E73", barcolor = "#009E73", ncp = 8) + 
  ylim(0,70) + 
  xlab("PCs") +
  my_theme +
  ggtitle("Explained variance Glasshouse")

which.max((get_pca_var(pca.Jena)$contrib)[,1]) # 1354
which.max((get_pca_var(pca.Jena)$contrib)[,2]) # 731
which.max((get_pca_var(pca.Jena)$contrib)[,3]) # 654
which.max((get_pca_var(pca.Jena)$contrib)[,4]) # 441

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

#SV_all_am = 100* apply(filter(R_3,Ctime=="am")[,51:2151],2,sd) / apply(filter(R_3,Ctime=="am")[,51:2151],2,mean)
#SV_all_am = t(as.data.frame(SV_all_am))
#colnames(SV_all_am) = 400:2500

#SV_all_noon = 100* apply(filter(R_3,Ctime=="noon")[,51:2151],2,sd) / apply(filter(R_3,Ctime=="noon")[,51:2151],2,mean)
#SV_all_noon = t(as.data.frame(SV_all_noon))
#colnames(SV_all_noon) = 400:2500

#SV_all_pm = 100* apply(filter(R_3,Ctime=="pm")[,51:2151],2,sd) / apply(filter(R_3,Ctime=="pm")[,51:2151],2,mean)
#SV_all_pm = t(as.data.frame(SV_all_pm))
#colnames(SV_all_pm) = 400:2500

# Absolute uncertainty
mean_AU_3 = t(apply(A_3[,51:2151],2,mean))
colnames(mean_AU_3) = 400:2500

# Relative uncertainty
RU_3 = 100 * A_3[,51:2151] / R_3[51:2151]
mean_RU_3 = apply(RU_3, 2, mean)
mean_RU_3 = t(as.data.frame(mean_RU_3))
colnames(mean_RU_3) = 400:2500

# Fig 5b
plot(as_spectra(SV_all_3), main="CV Field_UT",ylab = "CV (%)", xlab = "Wavelength (nm)",
     ylim = c(0,40), col = "#D55E00", lwd=3,xaxs = "i", yaxs = "i")
#areaplot(t(SV_all_3),col = rgb(252,78,7,max=255,alpha = 100),add=TRUE)
plot(as_spectra(SV_ref_3),col = "#000000", lwd=3,lty=2,add=TRUE)
plot(as_spectra(SV_TL_3),col = "#0072B2", lwd=3,lty=2,add=TRUE)
#plot(as_spectra(SV_all_am), col = "#fdbf6f",lwd=3,add=TRUE)
#plot(as_spectra(SV_all_noon), col = "#ff7f00",lwd=3,add=TRUE)
#plot(as_spectra(SV_all_pm), col = "#E69F00",lwd=3,add=TRUE)
plot(as_spectra(mean_RU_3), col = "#CC79A7",lwd=3,lty=3,add=TRUE)
abline(v=c(1000,1800), lty = 3)
plot_regions(as_spectra(SV_all_2),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add = TRUE)
legend(1900,15,legend=c("All","TLs","Ref","MU"),
       col=c("#D55E00","#0072B2","#000000","#CC79A7"), bg="transparent", 
       lty=c(1,2,2,3), lwd=3,cex=1,box.lty=0,seg.len = 2,x.intersp=0.5,y.intersp = 0.5)
Fig_5b <- recordPlot()

## PCA
R_4 = readRDS("~/Desktop/Nicotiana/Final/Data_rds/R_UT_dev.rds")
# remove colume Leaf_ID for now
R_4 = R_4[,-2154]
UT = rbind(R_3,R_4) %>%
        mutate(Genos = ifelse(Genotype_ID %in% c("UT_WT","pRESC2NC","pSOL3NC"), "Ref", "non-Ref"),
               GeneGroup = ifelse(Genos == "Ref","Ref","TL"))

pca.UT = prcomp(UT[,51:2151],scale. = TRUE)
# Fig S2e
Fig_S2e <- fviz_pca_ind(pca.UT, axes = c(1,2),label = "none", habillage = UT$Genos, palette = c("#D55E00", "#000000"),addEllipses = TRUE) + 
  my_theme +
  theme(legend.position = c(.11, .90))+
  ggtitle("PCA Field_UT")
# Fig 5f
Fig_5f <- fviz_pca_ind(pca.UT, axes = c(3,4),label = "none", habillage = UT$Genos, palette = c("#D55E00", "#000000"),addEllipses = TRUE) + 
  my_theme +
  theme(legend.position = c(.11, .90))+
  ggtitle("PCA Field_UT")
# fig S2b
Fig_S2b <- fviz_eig(pca.UT, addlabels = TRUE, main="Field_UT", barfill ="#D55E00", barcolor = "#D55E00", ncp = 8) + 
  ylim(0,75) + 
  xlab("PCs") +
  my_theme + 
  ggtitle("Explained variance Field_UT")

# PCA UT_lin
pca.UT.lin = prcomp(R_3[,51:2151],scale. = TRUE)

which.max((get_pca_var(pca.UT.lin)$contrib)[,1]) # 1516
which.max((get_pca_var(pca.UT.lin)$contrib)[,2]) # 1005
which.max((get_pca_var(pca.UT.lin)$contrib)[,3]) # 536
which.max((get_pca_var(pca.UT.lin)$contrib)[,4]) # 1801

# PCA UT_dev
pca.UT.dev = prcomp(R_4[,51:2151],scale. = TRUE)

which.max((get_pca_var(pca.UT.dev)$contrib)[,1]) # 1398
which.max((get_pca_var(pca.UT.dev)$contrib)[,2]) # 689
which.max((get_pca_var(pca.UT.dev)$contrib)[,3]) # 718
which.max((get_pca_var(pca.UT.dev)$contrib)[,4]) # 711

# Measurement of uncertainty comparision
# Figure 5d
plot(as_spectra(mean_RU_3), main="MU from (a)-(c)",col = "#D55E00",
     ylab = "Measurement of uncertainty (%)", xlab = "Wavelength",ylim = c(0,1.1),lwd=3, xaxs = "i", yaxs = "i")
plot(as_spectra(mean_RU),col = "#56B4E9", lwd=3,add=TRUE)
plot(as_spectra(mean_RU_2),col = "#009E73",lwd=3, add=TRUE)
abline(v=c(1000,1800), lty = 3)
plot_regions(as_spectra(mean_RU),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add = TRUE)
legend(1800,0.24,legend=c("Field_AZ","Field_UT","Glasshouse"),
       col=c("#56B4E9","#D55E00","#009E73"), bg="transparent",
       lty=c(1,1,1), lwd=c(3,3,3),cex=1,box.lty=0,seg.len = 2,x.intersp=0.5,y.intersp = 0.5)
Fig_5d <- recordPlot()

#absolute
#plot(as_spectra(mean_AU_3), main="Measurement uncertainty in three envs",col = "#D55E00",
#     ylab = "Measurement of uncertainty (%)", xlab = "Wavelength",ylim = c(0,0.003),lwd=3, xaxs = "i", yaxs = "i")
#plot(as_spectra(mean_AU),col = "#56B4E9", lwd=3,add=TRUE)
#plot(as_spectra(mean_AU_2),col = "#009E73",lwd=3, add=TRUE)
#abline(v=c(1000,1800), lty = 3)
#plot_regions(as_spectra(mean_RU),col = grDevices::rgb(0.7, 0.7, 0.7, 0.2), add = TRUE)

## Contributions to PC 1~4
spec_AZ = as_spectra(t(get_pca_var(pca.AZ)$contrib))
spec_Jena = as_spectra(t(get_pca_var(pca.Jena)$contrib))
spec_UT = as_spectra(t(get_pca_var(pca.UT)$contrib))

plot(spec_Jena[1,], 
     main="Contributions to PC1",
     xlab="Wavelength (nm)", ylab = "Contribution(%)",xaxs = "i", 
     col = "#009E73",lwd=3)
plot(spec_AZ[1,],lwd=3,col="#56B4E9",add=TRUE)
plot(spec_UT[1,],lwd=3,col="#D55E00",add=TRUE)
plot_regions(spec_Jena, col = grDevices::rgb(0.7, 0.7, 0.7, 0.2),add=TRUE)
abline(v=c(1000,1800), lty = 3)
legend(1800,0.02,legend=c("Field_AZ", "Field_UT","Glasshouse"),
       col=c("#56B4E9","#D55E00","#009E73"), bg="transparent", lty=1, lwd=3,cex=1,box.lty=0,seg.len = 2,x.intersp=0.5, y.intersp=0.5)
Fig_5h <- recordPlot()

plot(spec_AZ[2,], 
     main="Contributions to PC2",
     xlab="Wavelength (nm)", ylab = "Contribution(%)",xaxs = "i", 
     col = "#56B4E9",lwd=3)
plot(spec_Jena[2,],lwd=3,col="#009E73",add=TRUE)
plot(spec_UT[2,],lwd=3,col="#D55E00",add=TRUE)
plot_regions(spec_Jena, col = grDevices::rgb(0.7, 0.7, 0.7, 0.2),add=TRUE)
abline(v=c(1000,1800), lty = 3)
legend(1800,0.17,legend=c("Field_AZ", "Field_UT","Glasshouse"),
       col=c("#56B4E9","#D55E00","#009E73"), bg="transparent", lty=1, lwd=3,cex=1,box.lty=0,seg.len = 2,x.intersp=0.5, y.intersp=0.5)
Fig_5i <- recordPlot()

plot(spec_Jena[3,],col = "#009E73",
     main="Contributions to PC3",
     xlab="Wavelength (nm)", ylab = "Contribution(%)",xaxs = "i", 
     cex.lab=1.4, cex.axis=1.4,lwd=3, cex.main=1.5)
plot(spec_AZ[3,],col="#56B4E9",lwd=3,add=TRUE)
plot(spec_UT[3,],col="#D55E00",lwd=3,add=TRUE)
plot_regions(spec_Jena, col = grDevices::rgb(0.7, 0.7, 0.7, 0.2),add=TRUE)
abline(v=c(1000,1800), lty = 3)
legend(1800,0.3,legend=c("Field_AZ", "Field_UT","Glasshouse"),
       col=c("#56B4E9","#D55E00","#009E73"), bg="transparent", lty=1, lwd=3,cex=1,box.lty=0,seg.len = 2,x.intersp=0.5, y.intersp=0.5)
Fig_5j <- recordPlot()

plot(spec_AZ[4,],col="#56B4E9",
     main="Contributions to PC4",
     xlab="Wavelength (nm)", ylab = "Contribution(%)",xaxs = "i", 
     cex.lab=1.4, cex.axis=1.4,lwd=3, cex.main=1.5)
plot(spec_UT[4,],col="#D55E00",lwd=3,add=TRUE)
plot(spec_Jena[4,],col = "#009E73",lwd=3,add=TRUE)
plot_regions(spec_Jena, col = grDevices::rgb(0.7, 0.7, 0.7, 0.2),add=TRUE)
abline(v=c(1000,1800), lty = 3)
legend(1800,0.28,legend=c("Field_AZ", "Field_UT","Glasshouse"),
       col=c("#56B4E9","#D55E00","#009E73"), bg="transparent", lty=1, lwd=3,cex=1,box.lty=0,seg.len = 2,x.intersp=0.5, y.intersp=0.5)
Fig_5k <- recordPlot()

## All PCA
All = as.data.frame(rbind(AZ[,c(51:2152,2158:2160)],Jena[,c(51:2152,2156:2158)],UT[,c(51:2152,2155:2157)]))
pca.all = prcomp(All[,1:2101],scale. = TRUE)

Fig_3a <- fviz_pca_ind(pca.all, label = "none", title = "PCA across all samples", habillage = All$Site, palette = c("#56B4E9", "#D55E00","#009E73"),addEllipses = TRUE) + 
  my_theme +
  theme(legend.position = c(.13, .15))

Fig_3b <- fviz_pca_ind(pca.all, axes = c(3,4),label = "none", habillage = All$Site, palette = c("#56B4E9", "#D55E00","#009E73"),addEllipses = TRUE) + 
  my_theme +
  theme(legend.position = c(.13, .15))+
  ggtitle("")

# PCA UT_WT
All_UT = filter(All,Genotype_ID %in% c("UT-WT","UT_WT"))
pca.all.ut = prcomp(All_UT[,1:2101],scale. = TRUE)
Fig_3d <- fviz_pca_ind(pca.all.ut, label = "none", title = "PCA across UT-WT samples",habillage = All_UT$Site, palette = c("#56B4E9", "#D55E00","#009E73"), addEllipses = TRUE) + 
  my_theme +
  theme(legend.position = c(.13, .15))

Fig_3e <- fviz_pca_ind(pca.all.ut, axes = c(3,4),label = "none", habillage = All_UT$Site, palette = c("#56B4E9", "#D55E00","#009E73"), addEllipses = TRUE) + 
  my_theme +
  theme(legend.position = c(.13, .90))+
  ggtitle("")
# Contribution of variable to dimension 1~4, PCA_All
spec_all = as_spectra(t(get_pca_var(pca.all)$contrib))
plot(spec_all[4,], 
     main="Contributions to PCs 1-4",xaxs = "i",
     xlab="Wavelength (nm)", ylab = "Contribution(%)",
     col = "#000000",lwd=3)
plot(spec_all[1,],lwd=3,col="#56B4E9",add=TRUE)
plot(spec_all[2,],lwd=3,col="#D55E00",add=TRUE)
plot(spec_all[3,],lwd=3,col="#009E73",add=TRUE)
plot_regions(spec_all[1,], col = grDevices::rgb(0.7, 0.7, 0.7, 0.2),add=TRUE)
legend(2000,0.4,legend=c("PC 1", "PC 2","PC 3","PC 4"),
       col=c("#56B4E9","#D55E00","#009E73","#000000"), bg="transparent", lty=1, lwd=3,cex=1,box.lty=0,seg.len = 2,x.intersp=0.5, y.intersp=0.5)
Fig_3c <- recordPlot()

# Contribution of variable to dimension 1~4, PCA_UT-WT
spec_all_ut = as_spectra(t(get_pca_var(pca.all.ut)$contrib))
plot(spec_all_ut[4,], 
     main="Contributions to PCs 1-4",
     xlab="Wavelength (nm)", ylab = "Contribution(%)",xaxs = "i",
     col = "#000000",lwd=3)
plot(spec_all_ut[1,],lwd=3,col="#56B4E9",add=TRUE)
plot(spec_all_ut[2,],lwd=3,col="#D55E00",add=TRUE)
plot(spec_all_ut[3,],lwd=3,col="#009E73",add=TRUE)
plot_regions(spec_all_ut[1,], col = grDevices::rgb(0.7, 0.7, 0.7, 0.2),add=TRUE)
legend(2000,0.37,legend=c("PC 1", "PC 2","PC 3","PC 4"),
       col=c("#56B4E9","#D55E00","#009E73","#000000"), bg="transparent", lty=1, lwd=3,cex=1,box.lty=0,seg.len = 2,x.intersp=0.5, y.intersp=0.5)
Fig_3f <- recordPlot()

setwd("~/Desktop/Nicotiana/Final/Fig/")

# Fig 3
pdf("Fig3.pdf",width = 18,height = 12)
plot_grid(Fig_3a,Fig_3b,Fig_3c, Fig_3d, Fig_3e, Fig_3f, 
          nrow = 2,labels = c('(a)','(b)','(c)','(d)','(e)','(f)'))
dev.off()

# Fig 5
pdf("Fig5.pdf",width = 24,height = 18)
plot_grid(Fig_5a,Fig_5b,Fig_5c, Fig_5d, Fig_5e, Fig_5f,Fig_5g, NULL, Fig_5h, Fig_5i, Fig_5j,Fig_5k, 
          nrow = 3,labels = c('(a)','(b)','(c)','(d)','(e)','(f)','(g)','','(h)','(i)','(j)','(k)'))
dev.off()

# Fig S2
pdf("SFig2.pdf",width = 18,height = 12)
plot_grid(Fig_S2a,Fig_S2b,Fig_S2c, Fig_S2d, Fig_S2e, Fig_S2f, 
          nrow = 2,labels = c('(a)','(b)','(c)','(d)','(e)','(f)'))
dev.off()
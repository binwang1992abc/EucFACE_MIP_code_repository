

#############################################################################
#### Modelling CO2 x P at EucFACE: multi-model intercomparison           ####
#### Trace fate of the added P and its response to eCO2                  ####
#### Created by Bin Wang, 2025-04-27                                     ####
#### Email: wangbin1992@zju.edu.cn                                       ####
#############################################################################




### clear space
rm(list = ls());




#### preparations
#============================================
## packages
library(doBy)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(gridExtra)
library(patchwork)
library(reshape2)

## climate
climate.scenario <- 'FIX' # fixed climate in wet

## set pathways of input and output
setwd('D:/Analysis/EucFACE_MIP/EucFACE_MIP_P_x_CO2/')

## P gradients colours
Diverge_hsv_Palette <- colorspace::diverge_hcl(8)
Pall_color <- c( '#bbc258', '#FBEBBC',  '#75A4C9','#9D69B1')

## model names and labels
model.names <- c("A_GDAYP",
                 "B_ELMV1",
                 "C_CABLP",
                 "D_LPJGP",
                 "E_OCHDP",
                 "F_QUINC",
                 "G_OCHDX",
                 "H_QUJSM")
model.labels <- c("A_GDAYP" = "GDAYP",
                  "B_ELMV1" = "ELMV1",
                  "C_CABLP" = "CABLP",
                  "D_LPJGP" = "LPJGP",
                  "E_OCHDP" = "OCHDP",
                  "F_QUINC" = "QUINC",
                  "G_OCHDX" = "OCHDX",
                  "H_QUJSM" = "QUJSM")




### data preparation
# =======================================================
## read annual datasets
#-----------------------------------------------------------------------------------
ambNOP <- readRDS(paste0(getwd(),
                         "/data/compile_output/MIP_ALL_", 
                         climate.scenario, "_NOP_AMB_annual.rds"))
ambMDP <- readRDS(paste0(getwd(),
                         "/data/compile_output/MIP_ALL_", 
                         climate.scenario, "_MDP_AMB_annual.rds"))
ambHIP <- readRDS(paste0(getwd(),
                         "/data/compile_output/MIP_ALL_", 
                         climate.scenario, "_HIP_AMB_annual.rds"))

eleNOP <- readRDS(paste0(getwd(),
                         "/data/compile_output/MIP_ALL_", 
                         climate.scenario, "_NOP_ELE_annual.rds"))
eleMDP <- readRDS(paste0(getwd(),
                         "/data/compile_output/MIP_ALL_", 
                         climate.scenario, "_MDP_ELE_annual.rds"))
eleHIP <- readRDS(paste0(getwd(),
                         "/data/compile_output/MIP_ALL_", 
                         climate.scenario, "_HIP_ELE_annual.rds"))

## add P treatment, combine
ambNOP$PTRT <- eleNOP$PTRT <- "NOP"
ambMDP$PTRT <- eleMDP$PTRT <- "MDP"
ambHIP$PTRT <- eleHIP$PTRT <- "HIP"

ambDF <- rbind(ambNOP, rbind(ambMDP, ambHIP))
eleDF <- rbind(eleNOP, rbind(eleMDP, eleHIP))

## ignore two C-N models
ambDF <- ambDF[ambDF$ModName%in%c("A_GDAYP", "B_ELMV1", "C_CABLP", "D_LPJGP", "E_OCHDP",
                                  "F_QUINC", "G_OCHDX", "H_QUJSM"),]
eleDF <- eleDF[eleDF$ModName%in%c("A_GDAYP", "B_ELMV1", "C_CABLP", "D_LPJGP", "E_OCHDP",
                                  "F_QUINC", "G_OCHDX", "H_QUJSM"),]
## extract the P-related variables
fluxamb <- ambDF[,c("ModName", "PTRT", "YEAR", "PLITIN", "PCRLIN", "PFRLIN", "PWLIN",
                    "PUP", "PGMIN", "PMIN", "PBIOCHMIN", "PLEACH", "PGL", "PGW", "PGCR",
                    "PGFR", "PLRETR", "PWRETR", "PCRRETR", "PFRRETR", "PWEA", "PDEP",
                    "PFERT", 
                    "deltaPL", "deltaPW", "deltaPFR", "deltaPCR", "deltaPSTOR", "deltaPSOIL",
                    "deltaPPMIN", "deltaPLAB", "deltaPSEC", "deltaPOCC", "deltaPPAR",
                    "deltaPFLITA", "deltaPFLITB", "deltaPCLITB", "deltaPFLIT", "deltaPPORG")]

fluxele <- eleDF[,c("ModName", "PTRT", "YEAR", "PLITIN", "PCRLIN", "PFRLIN", "PWLIN",
                    "PUP", "PGMIN", "PMIN", "PBIOCHMIN", "PLEACH", "PGL", "PGW", "PGCR",
                    "PGFR", "PLRETR", "PWRETR", "PCRRETR", "PFRRETR", "PWEA", "PDEP",
                    "PFERT", 
                    "deltaPL", "deltaPW", "deltaPFR", "deltaPCR", "deltaPSTOR", "deltaPSOIL",
                    "deltaPPMIN", "deltaPLAB", "deltaPSEC", "deltaPOCC", "deltaPPAR",
                    "deltaPFLITA", "deltaPFLITB", "deltaPCLITB", "deltaPFLIT", "deltaPPORG")]

## prepare difference DF
fluxDF1 <- fluxamb[fluxamb$YEAR>=2020 & fluxamb$YEAR<=2022,]
fluxDF2 <- fluxele[fluxele$YEAR>=2020 & fluxele$YEAR<=2022,]

stDF1 <- fluxDF1[fluxDF1$YEAR==2020,]
stDF2 <- fluxDF2[fluxDF2$YEAR==2020,]

## sum over the three years
d2 <- dim(fluxDF1)[2]
for (i in unique(fluxDF1$ModName)) {
  for (j in unique(fluxDF1$PTRT)) {
    stDF1[stDF1$ModName==i & stDF1$PTRT==j,4:d2] <- colSums(fluxDF1[fluxDF1$ModName==i & fluxDF1$PTRT==j,4:d2], na.rm=T)
    stDF2[stDF2$ModName==i & stDF2$PTRT==j,4:d2] <- colSums(fluxDF2[fluxDF2$ModName==i & fluxDF2$PTRT==j,4:d2], na.rm=T)
  }
}




### data processing
# ============================================================================

## revise model-specific mass balance calculations
## because different models made different assumptions on which flux contributed to the overall P cycling

## PFERT in ELMV1 is doubled accounted for, hence making it zero
stDF1$PFERT[stDF1$ModName=="B_ELMV1"] <- 0.0
stDF2$PFERT[stDF2$ModName=="B_ELMV1"] <- 0.0

## additional variable PGOCL in LPJGP
stDF1$PGOCL <- 0.0
stDF2$PGOCL <- 0.0

stDF1$PGOCL[stDF1$ModName=="D_LPJGP"] <- with(stDF1[stDF1$ModName=="D_LPJGP",], PDEP+PWEA+PFERT-PLEACH-(deltaPL+deltaPW+deltaPCR+deltaPFR+deltaPSOIL+deltaPFLIT+deltaPCLITB))
stDF2$PGOCL[stDF2$ModName=="D_LPJGP"] <- with(stDF2[stDF2$ModName=="D_LPJGP",], PDEP+PWEA+PFERT-PLEACH-(deltaPL+deltaPW+deltaPCR+deltaPFR+deltaPSOIL+deltaPFLIT+deltaPCLITB))

## PFERT is doubled accounted for, hence making it zero
stDF1$PFERT[stDF1$ModName=="E_OCHDP"] <- 0.0
stDF2$PFERT[stDF2$ModName=="E_OCHDP"] <- 0.0

## P allocation to fruit and seed pool not reported, add them back
stDF1$deltaPOTHER <- 0.0
stDF2$deltaPOTHER <- 0.0

stDF1$deltaPOTHER[stDF1$ModName=="F_QUINC"] <- with(stDF1[stDF1$ModName=="F_QUINC",], PDEP+PWEA+PFERT-PLEACH-(deltaPL+deltaPW+deltaPCR+deltaPFR+deltaPSOIL+deltaPFLIT+deltaPCLITB))
stDF2$deltaPOTHER[stDF2$ModName=="F_QUINC"] <- with(stDF2[stDF2$ModName=="F_QUINC",], PDEP+PWEA+PFERT-PLEACH-(deltaPL+deltaPW+deltaPCR+deltaPFR+deltaPSOIL+deltaPFLIT+deltaPCLITB))

## PFERT is doubled accounted for, hence making it zero
stDF1$PFERT[stDF1$ModName=="G_OCHDX"] <- 0.0
stDF2$PFERT[stDF2$ModName=="G_OCHDX"] <- 0.0

## P allocation to fruit and seed pool not reported, add them back
stDF1$deltaPOTHER[stDF1$ModName=="H_QUJSM"] <- with(stDF1[stDF1$ModName=="H_QUJSM",], PDEP+PWEA+PFERT-PLEACH-(deltaPL+deltaPW+deltaPCR+deltaPFR+deltaPSOIL+deltaPFLIT+deltaPCLITB))
stDF2$deltaPOTHER[stDF2$ModName=="H_QUJSM"] <- with(stDF2[stDF2$ModName=="H_QUJSM",], PDEP+PWEA+PFERT-PLEACH-(deltaPL+deltaPW+deltaPCR+deltaPFR+deltaPSOIL+deltaPFLIT+deltaPCLITB))

## PFERT is a pool in GDAYP, need to revise it into a flux
## also need to account for the change in pool in PFERT (make it part of PSOIL)
stDF1$PFERT[stDF1$ModName=="A_GDAYP" & stDF1$PTRT=="MDP"] <- 4.5
stDF1$PFERT[stDF1$ModName=="A_GDAYP" & stDF1$PTRT=="HIP"] <- 9.0

stDF2$PFERT[stDF2$ModName=="A_GDAYP" & stDF2$PTRT=="MDP"] <- 4.5
stDF2$PFERT[stDF2$ModName=="A_GDAYP" & stDF2$PTRT=="HIP"] <- 9.0

stDF1$deltaPSOIL[stDF1$ModName=="A_GDAYP"] <- with(stDF1[stDF1$ModName=="A_GDAYP",], PDEP+PWEA+PFERT-PLEACH-(deltaPL+deltaPW+deltaPCR+deltaPFR+deltaPFLIT+deltaPCLITB))
stDF2$deltaPSOIL[stDF2$ModName=="A_GDAYP"] <- with(stDF2[stDF2$ModName=="A_GDAYP",], PDEP+PWEA+PFERT-PLEACH-(deltaPL+deltaPW+deltaPCR+deltaPFR+deltaPFLIT+deltaPCLITB))

## split soil P into labile and non-labile
stDF1$deltaPSOILOTHER <- with(stDF1, deltaPSOIL-deltaPLAB)
stDF2$deltaPSOILOTHER <- with(stDF2, deltaPSOIL-deltaPLAB)

## merge
stDF1$CO2TRT <- "AMB"
stDF2$CO2TRT <- "ELE"
myDF <- rbind(stDF1, stDF2)


## normalized P pools with P input
## ----------------------------------------------------------------------------------------------------------
## calculate P input, change in vegetation, out
myDF$PIN <- with(myDF, PDEP+PWEA+PFERT)
myDF$deltaPVEG <- with(myDF, deltaPL+deltaPW+deltaPCR+deltaPFR+deltaPOTHER+deltaPSTOR+deltaPFLIT+deltaPCLITB)
myDF$POUT <- with(myDF, PLEACH+PGOCL)

## normalize the 4 pools by their sum
myDF$P4Psum <- with(myDF, abs(POUT)+abs(deltaPVEG)+abs(deltaPLAB)+abs(deltaPSOILOTHER), na.rm=TRUE)
myDF$POUT_norm <- with(myDF, abs(POUT)/P4Psum)
myDF$deltaPVEG_norm <- with(myDF, abs(deltaPVEG)/P4Psum)
myDF$deltaPLAB_norm <- with(myDF, abs(deltaPLAB)/P4Psum)
myDF$deltaPSOILOTHER_norm <- with(myDF, abs(deltaPSOILOTHER)/P4Psum)

## splite P pools, fraction
subDF_4P <- myDF[,c("ModName", "PTRT", "YEAR", "CO2TRT", "POUT_norm", "deltaPVEG_norm", "deltaPLAB_norm",
                    "deltaPSOILOTHER_norm")]
plotDF_4P <- melt(subDF_4P, id.var=c("ModName", "PTRT", "YEAR", "CO2TRT"))
plotDF_4P$PTRT <- factor(plotDF_4P$PTRT, levels = c('NOP', 'MDP', 'HIP'))
plotDF_4P$variable <- factor(plotDF_4P$variable, levels = c("deltaPVEG_norm", "deltaPLAB_norm", "deltaPSOILOTHER_norm", "POUT_norm"))




### part 2, plot pattern, multi-individual models
# ===================================================

## plot fate of P (fraction) with P x CO2, 4 pools, by CO2TRT
plotDF_4P_1 <- subset(plotDF_4P, CO2TRT == 'AMB' & PTRT == 'HIP')
plotDF_4P_2 <- subset(plotDF_4P, CO2TRT == 'ELE' & PTRT == 'HIP')
plotDF_4P_3 <- plotDF_4P_1
plotDF_4P_3$value <- plotDF_4P_2$value - plotDF_4P_1$value 

# aCO2
p1 <-
  ggplot(data=plotDF_4P_1, 
         aes(x=ModName, y=value))+
  geom_bar(stat = "identity",
           aes(fill=variable),
           position=position_stack(),
           width = 0.8,
           col='grey30')+
  geom_hline(yintercept=c(0))+
  scale_fill_manual(name='Destination',
                    values=c("deltaPVEG_norm"=Pall_color[1],
                             "deltaPLAB_norm"=Pall_color[2],
                             "deltaPSOILOTHER_norm"=Pall_color[3],
                             "POUT_norm"=Pall_color[4]),
                    labels=c("deltaPVEG_norm"=expression(P[veg]),
                             "deltaPLAB_norm"=expression(P[labile]),
                             "deltaPSOILOTHER_norm"=expression(P[nonlabile]),
                             "POUT_norm"=expression(P[out])))+
  scale_x_discrete(name="",
                   # guide = guide_axis(n.dodge = 2),
                   labels=c("A_GDAYP" = "GDAYP",
                            "B_ELMV1" = "ELMV1",
                            "C_CABLP" = "CABLP",
                            "D_LPJGP" = "LPJGP",
                            "E_OCHDP" = "OCHDP",
                            "F_QUINC" = "QUINC",
                            "G_OCHDX" = "OCHDX",
                            "H_QUJSM" = "QUJSM"))+
  scale_y_continuous(limits = c(0, 1.1),
                     breaks = seq(0,1.0,0.2))+
  xlab(NULL)+
  ylab(expression("Allocation of extra P (fraction)"))+
  annotate('text', x=0.8, y=1.08, label='(a)', size=12)+
  # annotate('text', x=2, y=1.17, label='MDP', size=8)+
  theme_linedraw() +
  theme(panel.grid.minor=element_blank(),
        plot.margin = unit(c(0.5,1,1,0.8), 'cm'),
        axis.text.x=element_text(size=24, angle = 35, vjust = 1, hjust = 1),
        axis.title.x=element_text(size=24),
        axis.text.y=element_text(size=24),
        axis.title.y=element_text(size=24, vjust=2.5),
        legend.text=element_text(size=28),
        legend.title=element_text(size=28),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_blank(),
        legend.position="top",
        legend.box = 'horizontal',
        legend.box.just = 'left',
        legend.key.size = unit(2, 'cm'),
        legend.spacing.x = unit(0.5, 'cm'),
        strip.text = element_text(size=28),
        plot.title = element_text(size=24, face="bold"))

# CO2 effect
p2 <-
  ggplot(data=plotDF_4P_3, 
         aes(x=ModName, y=value))+
  geom_bar(stat = "identity",
           aes(fill=variable),
           position=position_stack(),
           width = 0.8,
           col='grey30')+
  geom_hline(yintercept=c(0))+
  scale_fill_manual(name='Destination',
                    values=c("deltaPVEG_norm"=Pall_color[1],
                             "deltaPLAB_norm"=Pall_color[2],
                             "deltaPSOILOTHER_norm"=Pall_color[3],
                             "POUT_norm"=Pall_color[4]),
                    labels=c("deltaPVEG_norm"=expression(P[veg]),
                             "deltaPLAB_norm"=expression(P[labile]),
                             "deltaPSOILOTHER_norm"=expression(P[nonlabile]),
                             "POUT_norm"=expression(P[out])))+
  scale_x_discrete(name="",
                   # guide = guide_axis(n.dodge = 2),
                   labels=c("A_GDAYP" = "GDAYP",
                            "B_ELMV1" = "ELMV1",
                            "C_CABLP" = "CABLP",
                            "D_LPJGP" = "LPJGP",
                            "E_OCHDP" = "OCHDP",
                            "F_QUINC" = "QUINC",
                            "G_OCHDX" = "OCHDX",
                            "H_QUJSM" = "QUJSM"))+
  # scale_y_continuous(limits = c(-0.22, 0.22),
  #                    breaks = seq(-0.2,0.2,0.1),
  #                    labels = c(-0.2, -0.1, 0, 0.10, 0.2))+
  xlab(NULL)+
  ylab(expression(CO[2] * " effect on extra P allocation (fraction)"))+
  annotate('text', x=0.8, y=0.22, label='(b)', size=12)+
  # annotate('text', x=2, y=1.17, label='MDP', size=8)+
  theme_linedraw() +
  theme(panel.grid.minor=element_blank(),
        plot.margin = unit(c(0.5,1,1,0.8), 'cm'),
        axis.text.x=element_text(size=24, angle = 35, vjust = 1, hjust = 1),
        axis.title.x=element_text(size=24),
        axis.text.y=element_text(size=24),
        axis.title.y=element_text(size=24, vjust=2.5),
        legend.text=element_text(size=28),
        legend.title=element_text(size=28),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_blank(),
        legend.position="top",
        legend.box = 'horizontal',
        legend.box.just = 'left',
        legend.key.size = unit(2, 'cm'),
        legend.spacing.x = unit(0.5, 'cm'),
        strip.text = element_text(size=28),
        plot.title = element_text(size=24, face="bold"))

## combine
(p1 + p2) +
  plot_layout(guides = 'collect') &
  theme(legend.position = 'top')



















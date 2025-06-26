

#############################################################################
#### Modelling CO2 x P at EucFACE: multi-model intercomparison           ####
#### Main and interactive effects on C and P processes                   ####
#### Created by Bin Wang, 2025-04-27                                     ####
#### Email: wangbin1992@zju.edu.cn                                       ####
#############################################################################


rm(list = ls());




#### basic preparation
# ============================================
## packages
library(doBy)
library(ggplot2)
library(ggpattern)
library(cowplot)
library(gridExtra)
library(patchwork)
library(reshape2)
library(RColorBrewer)
library(agricolae)

## pathway set
setwd('D:/Analysis/EucFACE_MIP/EucFACE_MIP_P_x_CO2/')

## color set
model.cols <- c("#ADAFB1", "#20a486ff", "#ECB884", "#4758A2", "#6CBEC3", "#619CD9",  "#E08D8B", "#AF8CBB")

## define climate
climate.scenario <- 'FIX' # fixed climate in wet

## define model names
model.labels <- c("A_GDAYP" = "GDAYP",
                  "B_ELMV1" = "ELMV1",
                  "C_CABLP" = "CABLP",
                  "D_LPJGP" = "LPJGP",
                  "E_OCHDP" = "OCHDP",
                  "F_QUINC" = "QUINC",
                  "G_OCHDX" = "OCHDX",
                  "H_QUJSM" = "QUJSM")



### data preparation
# ===================================================================================
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

## fixed variable, so ignore climate variables
ambDF[,c("PAR", "TAIR", "TSOIL", "VPD", "PREC")] <- NULL
eleDF[,c("PAR", "TAIR", "TSOIL", "VPD", "PREC")] <- NULL

## ignore NAs
ambDF[ambDF<=-999] <- NA
eleDF[eleDF<=-999] <- NA

## remove N-contained model
ambDF <- subset(ambDF, !(ModName == 'I_GDAYN' | ModName == 'J_LPJGN'))
eleDF <- subset(eleDF, !(ModName == 'I_GDAYN' | ModName == 'J_LPJGN'))

## remove 2012
ambDF <- subset(ambDF, YEAR >= 2013)
eleDF <- subset(eleDF, YEAR >= 2013)



## calculate variabale 
# -----------------------------------------------------------------------------
ambDF$CUE <- ambDF$NPP/ambDF$GPP # plant CUE
ambDF$Cveg <- rowSums(data.frame(ambDF$CL, # plant biomass
                                 ambDF$CW,
                                 ambDF$CFR,
                                 ambDF$CCR,
                                 ambDF$CSTOR),
                      na.rm = T)
ambDF$deltaCveg <- rowSums(data.frame(ambDF$deltaCL, # plant biomass growth
                                      ambDF$deltaCW,
                                      ambDF$deltaCFR,
                                      ambDF$deltaCCR,
                                      ambDF$deltaCSTOR),
                           na.rm = TRUE)
eleDF$CUE <- eleDF$NPP/eleDF$GPP
eleDF$Cveg <- rowSums(data.frame(eleDF$CL,
                                 eleDF$CW,
                                 eleDF$CFR,
                                 eleDF$CCR,
                                 eleDF$CSTOR),
                      na.rm = T)
eleDF$deltaCveg <- rowSums(data.frame(eleDF$deltaCL,
                                      eleDF$deltaCW,
                                      eleDF$deltaCFR,
                                      eleDF$deltaCCR,
                                      eleDF$deltaCSTOR),
                           na.rm = TRUE)

ambDF$CPL <- ambDF$CL/ambDF$PL # leaf C:P ratio
eleDF$CPL <- eleDF$CL/eleDF$PL
ambDF$CNL <- ambDF$CL/ambDF$NL # leaf C:N ratio
eleDF$CNL <- eleDF$CL/eleDF$NL
ambDF$CPW <- ambDF$CW/ambDF$PW # wood C:P ratio
eleDF$CPW <- eleDF$CW/eleDF$PW
ambDF$CPFR <- ambDF$CFR/ambDF$PFR # fine root C:P ratio
eleDF$CPFR <- eleDF$CFR/eleDF$PFR

ambDF[is.na(ambDF)] <- 0
eleDF[is.na(eleDF)] <- 0
ambDF$AG <- ambDF$CL + ambDF$CW + ambDF$CSTOR # aboveboard biomass
ambDF$BG <- ambDF$CCR + ambDF$CFR # belowground biomass
ambDF$ABR <- ambDF$AG/ambDF$BG  # above- to below-ground biomass allocation
eleDF$AG <- eleDF$CL + eleDF$CW + ele_DF$CSTOR
eleDF$BG <- eleDF$CCR + eleDF$CFR
eleDF$ABR <- eleDF$AG/eleDF$BG
ambDF$RSOIL <- ambDF$RCR + ambDF$RFR + ambDF$RHET  # soil respiration
eleDF$RSOIL <- eleDF$RCR + eleDF$RFR + eleDF$RHET


## organize to make neat
ambDF$CO2TRT <- 'ambient'
eleDF$CO2TRT <- 'elevated'
ambDF <- ambDF[c(1,2,163,151,3:150,152:162)]
eleDF <- eleDF[c(1,2,163,151,3:150,152:162)]




#### data processing
# ===================================================================

NOP_amb <- subset(ambDF, (YEAR >= 2020 & YEAR <= 2029) &
                    PTRT == 'NOP')
MDP_amb <- subset(ambDF, (YEAR >= 2020 & YEAR <= 2029) &
                    PTRT == 'MDP')
HIP_amb <- subset(ambDF, (YEAR >= 2020 & YEAR <= 2029) &
                    PTRT == 'HIP')
NOP_ele <- subset(eleDF, (YEAR >= 2020 & YEAR <= 2029) &
                    PTRT == 'NOP')
MDP_ele <- subset(eleDF, (YEAR >= 2020 & YEAR <= 2029) &
                    PTRT == 'MDP')
HIP_ele <- subset(eleDF, (YEAR >= 2020 & YEAR <= 2029) &
                    PTRT == 'HIP')


## calculate main and interactive effects of P and CO2
# ---------------------------------------------------------
d1 <- dim(NOP_amb)[2]
# magnitude
MDP_eff <- NOP_amb # MDP effect
MDP_eff[, 5:d1] <- MDP_amb[, 5:d1] - NOP_amb[, 5:d1]
MDP_eff$Effect <- 'MDP'
  
HIP_eff <- NOP_amb # HIP effect
HIP_eff[, 5:d1] <- HIP_amb[, 5:d1] - NOP_amb[, 5:d1]
HIP_eff$Effect <- 'HIP'

ele_eff <- NOP_amb # CO2 effect, NOP
ele_eff[, 5:d1] <- NOP_ele[, 5:d1] - NOP_amb[, 5:d1]
ele_eff$Effect <- 'CO2'

ele_MPeff <- MDP_amb # CO2 effect, MDP
ele_MPeff[, 5:d1] <- MDP_ele[, 5:d1] - MDP_amb[, 5:d1]
ele_MPeff$Effect <- 'CO2_MDP'

ele_HPeff <- HIP_amb # CO2 effect, HIP
ele_HPeff[, 5:d1] <- HIP_ele[, 5:d1] - HIP_amb[, 5:d1]
ele_HPeff$Effect <- 'CO2_HIP'

MPE_eff <- NOP_amb # MDP + CO2 effect
MPE_eff[, 5:d1] <- MDP_ele[, 5:d1] - NOP_amb[, 5:d1]
MPE_eff$Effect <- 'MDP+CO2'

HPE_eff <- NOP_amb # HIP + CO2 effect
HPE_eff[, 5:d1] <- HIP_ele[, 5:d1] - NOP_amb[, 5:d1]
HPE_eff$Effect <- 'HIP+CO2'

MPI_eff <- NOP_amb # MDP x CO2 effect
MPI_eff[, 5:d1] <- (MDP_ele[, 5:d1] - MDP_amb[, 5:d1]) - (NOP_ele[, 5:d1] - NOP_amb[, 5:d1])
MPI_eff$Effect <- 'MDP_x_CO2'

HPI_eff <- NOP_amb # HIP x CO2 effect
HPI_eff[, 5:d1] <- (HIP_ele[, 5:d1] - HIP_amb[, 5:d1]) - (NOP_ele[, 5:d1] - NOP_amb[, 5:d1])
HPI_eff$Effect <- 'HIP_x_CO2'


# percentage 
MDP_eff1 <- NOP_amb # MDP effect
MDP_eff1[, 5:d1] <- (MDP_amb[, 5:d1] - NOP_amb[, 5:d1])/NOP_amb[, 5:d1]*100
MDP_eff1$Effect <- 'MDP'

HIP_eff1 <- NOP_amb # HIP effect
HIP_eff1[, 5:d1] <- (HIP_amb[, 5:d1] - NOP_amb[, 5:d1])/NOP_amb[, 5:d1]*100
HIP_eff1$Effect <- 'HIP'

ele_eff1 <- NOP_amb # CO2 effect
ele_eff1[, 5:d1] <- (NOP_ele[, 5:d1] - NOP_amb[, 5:d1])/NOP_amb[, 5:d1]*100
ele_eff1$Effect <- 'CO2'

ele_MPeff1 <- MDP_amb # CO2 effect, MDP
ele_MPeff1[, 5:d1] <- (MDP_ele[, 5:d1] - MDP_amb[, 5:d1])/MDP_amb[, 5:d1]*100
ele_MPeff1$Effect <- 'CO2_MDP'

ele_HPeff1 <- HIP_amb # CO2 effect, HIP
ele_HPeff1[, 5:d1] <- (HIP_ele[, 5:d1] - HIP_amb[, 5:d1])/HIP_amb[, 5:d1]*100
ele_HPeff1$Effect <- 'CO2_HIP'

MPE_eff1 <- NOP_amb # MDP + CO2 effect
MPE_eff1[, 5:d1] <- (MDP_ele[, 5:d1] - NOP_amb[, 5:d1])/NOP_amb[, 5:d1]*100
MPE_eff1$Effect <- 'MDP+CO2'

HPE_eff1 <- NOP_amb # HIP + CO2 effect
HPE_eff1[, 5:d1] <- (HIP_ele[, 5:d1] - NOP_amb[, 5:d1])/NOP_amb[, 5:d1]*100
HPE_eff1$Effect <- 'HIP+CO2'

MPI_eff1 <- NOP_amb # MDP x CO2 effect
MPI_eff1[, 5:d1] <- (MDP_ele[, 5:d1] - MDP_amb[, 5:d1])/MDP_amb[, 5:d1]*100 - (NOP_ele[, 5:d1] - NOP_amb[, 5:d1])/NOP_amb[, 5:d1]*100
MPI_eff1$Effect <- 'MDP_x_CO2'

HPI_eff1 <- NOP_amb # HIP x CO2 effect
HPI_eff1[, 5:d1] <- (HIP_ele[, 5:d1] - HIP_amb[, 5:d1])/HIP_amb[, 5:d1]*100 - (NOP_ele[, 5:d1] - NOP_amb[, 5:d1])/NOP_amb[, 5:d1]*100
HPI_eff1$Effect <- 'HIP_x_CO2'


## combine data
PC_effect_mag <- rbind(MDP_eff, 
                       rbind(HIP_eff, 
                             rbind(ele_eff, 
                                   rbind(ele_MPeff, 
                                         rbind(ele_HPeff,
                                               rbind(MPE_eff, 
                                                     rbind(HPE_eff, 
                                                           rbind(MPI_eff,
                                                                 HPI_eff))))))))
PC_effect_per <- rbind(MDP_eff1, 
                       rbind(HIP_eff1, 
                             rbind(ele_eff1, 
                                   rbind(ele_MPeff1, 
                                         rbind(ele_HPeff1,
                                               rbind(MPE_eff1, 
                                                     rbind(HPE_eff1, 
                                                           rbind(MPI_eff1,
                                                                 HPI_eff1))))))))
PC_effect_mag$Value <- 'magnitude'
PC_effect_per$Value <- 'percentage'
PC_effect <- rbind(PC_effect_mag, PC_effect_per)
# add time periods, 2020-2022, 2023-2029
PC_effect$Period[PC_effect$YEAR >= 2020 & PC_effect$YEAR <= 2022] <- 'temp'
PC_effect$Period[PC_effect$YEAR >= 2023 & PC_effect$YEAR <= 2029] <- 'futu'

PE_effect <- subset(PC_effect, select = c('ModName', 'YEAR', 'Effect', 'Value', 'Period',
                                          'GPP', 'deltaCveg', 'NEP', 
                                          'CUE', 'NPP', 'RAU', 'ABR',
                                          'CPL', 'CPW', 'CPFR',  
                                          'RHET', 'RSOIL', 'CSOIL',
                                          'LAI', 'Cveg', 'CNL',
                                          'PUP', 'NUP', 'PLAB', 'PMIN'))

## Model-specific effects
PE_effect_model <- summaryBy(.~ModName+Effect+Value+Period, data = PE_effect[-2], FUN=c(mean, sd), keep.names = T, na.rm=TRUE)

## calculate multi-model mean
PE_effect_mip <- summaryBy(.~Effect+Value+Period, data = PE_effect_model[2:24], FUN=c(mean, sd), keep.names = T, na.rm=T)
names(PE_effect_mip)[4:43] <- c('GPP', 'deltaCveg', 'NEP', 
                                'CUE', 'NPP', 'RAU', 'ABR',
                                'CPL', 'CPW', 'CPFR',  
                                'RHET', 'RSOIL', 'CSOIL',
                                'LAI', 'Cveg', 'CNL',
                                'PUP', 'NUP', 'PLAB', 'PMIN',
                                paste(c('GPP', 'deltaCveg', 'NEP', 
                                        'CUE', 'NPP', 'RAU', 'ABR',
                                        'CPL', 'CPW', 'CPFR',  
                                        'RHET', 'RSOIL', 'CSOIL',
                                        'LAI', 'Cveg', 'CNL',
                                        'PUP', 'NUP', 'PLAB', 'PMIN'), 'sd', sep = '_'))
PE_effect_mip$ModName <- 'M_M'


## incorporate model-specific scores to effects based on Jiang et al. 2024 (Science Advances) 
Model_score <- data.frame(ModName = unique(PE_effect_model$ModName),
                          Score = c(13,9,8,11,14,14,9,16)/94)

PE_effect_model1 <- merge(PE_effect_model[(1:24)], Model_score, by = 'ModName', all.x = TRUE)
PE_effect_model1[5:24] <- PE_effect_model1[5:24]*PE_effect_model1$Score
# calcualte multi-model mean weighted by model score
PE_effect_mip1 <- summaryBy(.~Effect+Value+Period, data = PE_effect_model1[2:24], FUN=sum, keep.names = T, na.rm=F)
PE_effect_mip1_sd <- PE_effect_mip[c(1:3,24:43)]
PE_effect_mip1_sd[4:23] <- NA
PE_effect_mip1 <- merge(PE_effect_mip1, PE_effect_mip1_sd, by= c('Effect', 'Value', 'Period'))
PE_effect_mip1$ModName <- 'M_M1'
# arrange
PE_effect_mip <- PE_effect_mip[c(44,1:3,4:43)]
PE_effect_mip1 <- PE_effect_mip1[c(44,1:3,4:43)]
names(PE_effect_mip1) <- names(PE_effect_model) <- names(PE_effect_mip)
# combine all
PE_effect_all <- rbind(PE_effect_mip, rbind(PE_effect_mip1, PE_effect_model))
PE_effect_all$ModName <- factor(PE_effect_all$ModName, levels = c('M_M',
                                                                  'M_M1',
                                                                  'A_GDAYP',
                                                                  'B_ELMV1',
                                                                  'C_CABLP',
                                                                  'D_LPJGP',
                                                                  'E_OCHDP',
                                                                  'F_QUINC',
                                                                  'G_OCHDX',
                                                                  'H_QUJSM'))




## plot main and interactive effects of CO2 and HP, 2020-2022
# ==============================================================================================
PE_effect_plot <- subset(PE_effect_all, Period == 'temp' & !(Effect == 'HIP+CO2' |
                                                               Effect == 'CO2_MDP' |
                                                               Effect == 'MDP' | 
                                                               Effect == 'MDP+CO2' | 
                                                               Effect == 'MDP_x_CO2'))

## plot main and interactive effects of CO2 and HP
# -----------------------------------------------------------------------------------
{
  # GPP
  summary(aov(GPP~Effect, data = subset(PE_effect_plot, Value == 'percentage' &
                                          (Effect == 'CO2' | Effect == 'CO2_HIP'))))
  t.test(GPP~Effect, data = subset(PE_effect_plot, Value == 'percentage' &
                                     (Effect == 'CO2' | Effect == 'CO2_HIP')))
  pc1 <-
    ggplot(data = subset(PE_effect_plot, Value == 'percentage' &
                           (Effect == 'CO2' | Effect == 'CO2_HIP')), 
           aes(x=Effect, y=GPP))+
    stat_boxplot(geom = 'errorbar', position = position_dodge(width = 0.5), width = 0.5, size=1)+
    geom_boxplot(position = position_dodge(width = 0.5),
                 width =0.5,
                 linewidth=1,
                 outlier.shape = NA)+
    geom_point(aes(fill=ModName),position = position_jitter(0.2), size=7.5, pch=21, alpha=0.8)+
    scale_fill_manual(name="",
                      values=c('black', 'red', model.cols))+
    scale_x_discrete(name='',
                     labels=c('CO2' = 'aP',
                              'CO2_HIP' = 'eP'))+
    scale_y_continuous(limits = c(-10,33),
                       breaks = seq(-10,30,10))+
    geom_hline(yintercept = 0, lty=2)+
    ylab(expression(CO[2] * " effect on GPP (%)"  ))+
    annotate('text', x=0.65, y=31, label='(a)', size=10)+
    theme_linedraw() +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          plot.margin = unit(c(0.8,1.3,0,0.5), 'cm'),
          plot.title = element_blank(),
          axis.text.x=element_text(size=22, vjust = 0),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=22),
          axis.title.y=element_text(size=22, vjust = -0.3),
          axis.ticks.length = unit(0.2, 'cm'),
          legend.text=element_text(size=23),
          legend.title=element_text(size=20),
          legend.position='none',
          legend.box = 'horizontal',
          legend.box.just = 'left',
          legend.key.size = unit(1.8, 'cm'))
  
  # biomass growth
  summary(aov(deltaCveg~Effect, data = subset(PE_effect_plot, Value == 'magnitude' &
                                                (Effect == 'CO2' | Effect == 'CO2_HIP'))))
  t.test(deltaCveg~Effect, data = subset(PE_effect_plot, Value == 'magnitude' &
                                           (Effect == 'CO2' | Effect == 'CO2_HIP')))
  pc2 <-
    ggplot(data = subset(PE_effect_plot, Value == 'magnitude' &
                           (Effect == 'CO2' | Effect == 'CO2_HIP')), 
           aes(x=Effect, y=deltaCveg/1000))+
    stat_boxplot(geom = 'errorbar', position = position_dodge(width = 0.5), width = 0.5, size=1)+
    geom_boxplot(position = position_dodge(width = 0.5),
                 width =0.5,
                 linewidth=1,
                 outlier.shape = NA)+
    geom_point(aes(fill=ModName),position = position_jitter(0.2), size=7.5, pch=21, alpha=0.8)+
    scale_fill_manual(name="",
                      values=c('black', 'red', model.cols))+
    scale_x_discrete(name='',
                     labels=c('CO2' = 'aP',
                              'CO2_HIP' = 'eP'))+
    scale_y_continuous(limits = c(-0.2,0.3),
                       breaks = seq(-0.2,0.3,0.1))+
    geom_hline(yintercept = 0, lty=2)+
    ylab(expression(CO[2] * " effect on " * Delta * C[veg] * " (kg C " * m^-2 * " " * yr^-1 * ")"))+
    annotate('text', x=0.65, y=0.28, label='(b)', size=10)+
    theme_linedraw() +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          plot.margin = unit(c(0.8,1.3,0,0.5), 'cm'),
          plot.title = element_blank(),
          axis.text.x=element_text(size=22, vjust = 0),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=22),
          axis.title.y=element_text(size=22),
          axis.ticks.length = unit(0.2, 'cm'),
          legend.text=element_text(size=23),
          legend.title=element_text(size=20),
          legend.position='none',
          legend.box = 'horizontal',
          legend.box.just = 'left',
          legend.key.size = unit(1.8, 'cm'))
  
  # NEP
  summary(aov(NEP~Effect, data = subset(PE_effect_plot, Value == 'magnitude' &
                                          (Effect == 'CO2' | Effect == 'CO2_HIP'))))
  t.test(NEP~Effect, data = subset(PE_effect_plot, Value == 'magnitude' &
                                     (Effect == 'CO2' | Effect == 'CO2_HIP')))
  pc3 <-
    ggplot(data = subset(PE_effect_plot, Value == 'magnitude' &
                           (Effect == 'CO2' | Effect == 'CO2_HIP')), 
           aes(x=Effect, y=NEP/1000))+
    stat_boxplot(geom = 'errorbar', position = position_dodge(width = 0.5), width = 0.5, size=1)+
    geom_boxplot(position = position_dodge(width = 0.5),
                 width =0.5,
                 size=1,
                 outlier.shape = NA)+
    geom_point(aes(fill=ModName),position = position_jitter(0.2), size=7.5, pch=21, alpha=0.8)+
    scale_fill_manual(name="",
                      values=c('black', 'red', model.cols))+
    scale_x_discrete(name='',
                     labels=c('CO2' = 'aP',
                              'CO2_HIP' = 'eP'))+
    scale_y_continuous(limits = c(-0.2,0.3),
                       breaks = seq(-0.2,0.3,0.1))+
    geom_hline(yintercept = 0, lty=2)+
    ylab(expression(CO[2] * " effect on " * NEP * " (kg C " * m^-2 * " " * yr^-1 * ")"))+
    annotate('text', x=0.65, y=0.28, label='(c)', size=10)+
    theme_linedraw() +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          plot.margin = unit(c(0.8,1.3,0,0.5), 'cm'),
          plot.title = element_blank(),
          axis.text.x=element_text(size=22, vjust = 0),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=22),
          axis.title.y=element_text(size=22),
          axis.ticks.length = unit(0.2, 'cm'),
          legend.text=element_text(size=23),
          legend.title=element_text(size=20),
          legend.position='none',
          legend.box = 'horizontal',
          legend.box.just = 'left',
          legend.key.size = unit(1.8, 'cm'))
  
  
  
  ## interactive effect of CO2 and HP enrichment
  # ------------------------------------------------------------------------------------------
  # GPP
  pi1 <-
    ggplot(data = subset(PE_effect_plot, Effect == 'HIP_x_CO2' & Value == 'percentage'))+
    geom_bar(stat = "identity",
             mapping=aes(ModName, GPP, fill=ModName),
             position=position_dodge(preserve = 'single'),
             col="black",
             width = 0.6) +
    geom_errorbar(mapping=aes(x = ModName, 
                              ymax = GPP + GPP_sd,
                              ymin = GPP - GPP_sd,
                              col= 'black'),
                  position=position_dodge(preserve = 'single'),
                  col="black",
                  width = 0.4)+
    geom_vline(xintercept=c(2.5), lty=3, linewidth=1)+
    scale_fill_manual(name="",
                      values=c('grey30', 'red', model.cols),
                      labels=c('M_M', 'M_M1', model.labels))+
    scale_x_discrete(name="",
                     # guide = guide_axis(n.dodge = 2),
                     labels=c("A_GDAYP" = "GDAYP",
                              "B_ELMV1" = "ELMV1",
                              "C_CABLP" = "CABLP",
                              "D_LPJGP" = "LPJGP",
                              "E_OCHDP" = "OCHDP",
                              "F_QUINC" = "QUINC",
                              "G_OCHDX" = "OCHDX",
                              "H_QUJSM" = "QUJSM",
                              "M_M" = "M_M",
                              "M_M1" = 'W_M'))+
    scale_y_continuous(limits = c(-20, 10),
                       breaks = seq(-20,10,10))+
    geom_hline(yintercept=0, lty=1)+
    ylab(NULL)+
    ylab(expression("Interactive effect on " * GPP * " (%)"))+
    annotate('text', x=1, y=9, label='(d)', size=10)+
    theme_linedraw() +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          plot.margin = unit(c(1,1,0.2,0.5), 'cm'),
          plot.title = element_text(size=24, face="bold",
                                    hjust = 0.5),
          axis.text.x=element_text(size=20, angle = 38, vjust = 1.0, hjust = 1.0),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=22),
          axis.title.y=element_text(size=22, vjust = 1),
          axis.ticks.length = unit(0.2, 'cm'),
          legend.text=element_text(size=23),
          legend.title=element_text(size=20),
          legend.position='none',
          legend.box = 'horizontal',
          legend.box.just = 'left',
          legend.key.size = unit(1.3, 'cm'))
  
  # biomass growth
  pi2 <-
    ggplot(data = subset(PE_effect_plot, Effect == 'HIP_x_CO2' & Value == 'magnitude'))+
    geom_bar(stat = "identity",
             mapping=aes(ModName, deltaCveg/1000, fill=ModName),
             position=position_dodge(preserve = 'single'),
             col="black",
             width = 0.6) +
    geom_errorbar(mapping=aes(x = ModName, 
                              ymax = deltaCveg/1000 + deltaCveg_sd/1000,
                              ymin = deltaCveg/1000 - deltaCveg_sd/1000,
                              col= 'black'),
                  position=position_dodge(preserve = 'single'),
                  col="black",
                  width = 0.4)+
    geom_vline(xintercept=c(2.5), lty=3, linewidth=1)+
    scale_fill_manual(name="",
                      values=c('grey30', 'red', model.cols),
                      labels=c('M_M', 'M_M1', model.labels))+
    scale_x_discrete(name="",
                     # guide = guide_axis(n.dodge = 2),
                     labels=c("A_GDAYP" = "GDAYP",
                              "B_ELMV1" = "ELMV1",
                              "C_CABLP" = "CABLP",
                              "D_LPJGP" = "LPJGP",
                              "E_OCHDP" = "OCHDP",
                              "F_QUINC" = "QUINC",
                              "G_OCHDX" = "OCHDX",
                              "H_QUJSM" = "QUJSM",
                              "M_M" = "M_M",
                              "M_M1" = 'W_M'))+
    scale_y_continuous(limits = c(-0.61, 0.31),
                       breaks = seq(-0.6,0.2,0.2),
                       labels = c(-0.6,-0.4,-0.2,0,0.2))+
    geom_hline(yintercept=0, lty=1)+
    ylab(NULL)+
    ylab(expression("Interactive effect on " * Delta * C[veg] * " (kg C " * m^-2 * " " * yr^-1 * ")"))+
    annotate('text', x=1, y=0.28, label='(e)', size=10)+
    theme_linedraw() +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          plot.margin = unit(c(1,1,0.2,0.5), 'cm'),
          plot.title = element_text(size=24, face="bold",
                                    hjust = 0.5),
          axis.text.x=element_text(size=20, angle = 38, vjust = 1.0, hjust = 1.0),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=22),
          axis.title.y=element_text(size=22, vjust = 1),
          axis.ticks.length = unit(0.2, 'cm'),
          legend.text=element_text(size=23),
          legend.title=element_text(size=20),
          legend.position='none',
          legend.box = 'horizontal',
          legend.box.just = 'left',
          legend.key.size = unit(1.3, 'cm'))
  
  # NEP
  pi3 <-
    ggplot(data = subset(PE_effect_plot, Effect == 'HIP_x_CO2' & Value == 'magnitude'))+
    geom_bar(stat = "identity",
             mapping=aes(ModName, NEP/1000, fill=ModName),
             position=position_dodge(preserve = 'single'),
             col="black",
             width = 0.6) +
    geom_errorbar(mapping=aes(x = ModName, 
                              ymax = NEP/1000 + NEP_sd/1000,
                              ymin = NEP/1000 - NEP_sd/1000,
                              col= 'black'),
                  position=position_dodge(preserve = 'single'),
                  col="black",
                  width = 0.4)+
    geom_vline(xintercept=c(2.5), lty=3, linewidth=1)+
    scale_fill_manual(name="",
                      values=c('grey30', 'red', model.cols),
                      labels=c('M_M', 'M_M1', model.labels))+
    scale_x_discrete(name="",
                     # guide = guide_axis(n.dodge = 2),
                     labels=c("A_GDAYP" = "GDAYP",
                              "B_ELMV1" = "ELMV1",
                              "C_CABLP" = "CABLP",
                              "D_LPJGP" = "LPJGP",
                              "E_OCHDP" = "OCHDP",
                              "F_QUINC" = "QUINC",
                              "G_OCHDX" = "OCHDX",
                              "H_QUJSM" = "QUJSM",
                              "M_M" = "M_M",
                              "M_M1" = 'W_M'))+
    scale_y_continuous(limits = c(-0.40, 0.40),
                       breaks = seq(-0.4,0.4,0.2))+
    geom_hline(yintercept=0, lty=1)+
    ylab(NULL)+
    ylab(expression("Interactive effect on " * NEP * " (kg C " * m^-2 * " " * yr^-1 * ")"))+
    annotate('text', x=1, y=0.37, label='(f)', size=10)+
    theme_linedraw() +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          plot.margin = unit(c(1,1,0.2,0.5), 'cm'),
          plot.title = element_text(size=24, face="bold",
                                    hjust = 0.5),
          axis.text.x=element_text(size=20, angle = 38, vjust = 1.0, hjust = 1.0),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=22),
          axis.title.y=element_text(size=22, vjust = 1),
          axis.ticks.length = unit(0.2, 'cm'),
          legend.text=element_text(size=23),
          legend.title=element_text(size=20),
          legend.position='none',
          legend.box = 'horizontal',
          legend.box.just = 'left',
          legend.key.size = unit(1.3, 'cm'))
  
  
  ## relationship between predicted CO2 and P effects
  # ------------------------------------------------------------------------------------------
  P_CO2_plot1 <- subset(PE_effect_plot, Effect == 'CO2')
  P_CO2_plot2 <- subset(PE_effect_plot, Effect == 'HIP')
  names(P_CO2_plot1)[5:44] <- paste(names(P_CO2_plot1)[5:44], 'CO2', sep = '_')
  names(P_CO2_plot2)[5:44] <- paste(names(P_CO2_plot2)[5:44], 'P', sep = '_')
  P_CO2_plot <- cbind(P_CO2_plot1, P_CO2_plot2)
  
  # GPP
  summary(lm(GPP_CO2~GPP_P, data = subset(P_CO2_plot, Value == 'percentage')))
  p1 <-
    ggplot(data = subset(P_CO2_plot, Value == 'percentage'),
           aes(GPP_P, GPP_CO2))+
    geom_point(aes(fill=ModName), pch=21, size=8)+
    geom_errorbar(aes(x=GPP_P,
                      ymax = GPP_CO2 + GPP_sd_CO2,
                      ymin = GPP_CO2 - GPP_sd_CO2),
                  linewidth = 1.2,
                  col='grey40',
                  alpha=0.5,
                  width=1.5)+
    geom_errorbar(aes(y=GPP_CO2,
                      xmax = GPP_P + GPP_sd_P,
                      xmin = GPP_P - GPP_sd_P),
                  linewidth = 1.2,
                  col='grey40',
                  alpha=0.5,
                  width=1.5)+
    geom_smooth(method = 'lm', lty=2, col='blue', linewidth=2.5, se=T, alpha=0.2)+
    scale_x_continuous(limits = c(-20,50),
                       breaks = seq(-20,50,10))+
    scale_y_continuous(limits = c(0,40),
                       breaks = seq(0,40,10))+
    scale_fill_manual(name="",
                      values=c('grey30', 'red', model.cols),
                      labels=c('M_M', 'M_M1', model.labels))+
    geom_hline(yintercept = 0, lty=2)+
    geom_vline(xintercept = 0, lty=2)+
    ylab(expression(CO[2] * '  effect on GPP (%)'))+
    xlab('P effect on GPP (%)')+
    # ylab(expression("P effect (kg C " * m^-2 * " " * yr^-1 * ")"  ))+
    annotate('text', x=-18, y=38, label='(g)', size=10)+
    annotate('text', x=28, y=37, label=substitute(~R^2~"= 0.07,"~italic(p)~"= 0.47"), size=9, family = 'serif')+
    theme_linedraw() +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          plot.margin = unit(c(0.8,1,0,0.5), 'cm'),
          axis.ticks.length = unit(0.2, 'cm'),
          plot.title = element_blank(),
          axis.text.x=element_text(size=22, vjust = 0),
          axis.title.x=element_text(size=22, vjust = 4),
          # axis.title.x=element_blank(),
          axis.text.y=element_text(size=22),
          axis.title.y=element_text(size=22, vjust = -0.3),
          legend.text=element_text(size=22),
          legend.title=element_blank(),
          legend.position='none',
          legend.box = 'horizontal',
          legend.box.just = 'left',
          legend.key.size = unit(1.8, 'cm'))
  
  # biomass growth
  summary(lm(deltaCveg_CO2~deltaCveg_P, data = subset(P_CO2_plot, Value == 'magnitude')))
  p2 <-
    ggplot(data = subset(P_CO2_plot, Value == 'magnitude'),
           aes(deltaCveg_P, deltaCveg_CO2))+
    geom_point(aes(fill=ModName), pch=21, size=8)+
    geom_errorbar(aes(x=deltaCveg_P,
                      ymax = deltaCveg_CO2 + deltaCveg_sd_CO2,
                      ymin = deltaCveg_CO2 - deltaCveg_sd_CO2),
                  linewidth = 1.2,
                  col='grey40',
                  alpha=0.5,
                  width=1.5)+
    geom_errorbar(aes(y=deltaCveg_CO2,
                      xmax = deltaCveg_P + deltaCveg_sd_P,
                      xmin = deltaCveg_P - deltaCveg_sd_P),
                  linewidth = 1.2,
                  col='grey40',
                  alpha=0.5,
                  width=1.5)+
    geom_smooth(method = 'lm', lty=1, col='blue', linewidth=2.5, se=T, alpha=0.2)+
    scale_x_continuous(limits = c(-210,610),
                       breaks = seq(-200,600,200))+
    scale_y_continuous(limits = c(-140,340),
                       breaks = seq(-100,300,100))+
    scale_fill_manual(name="",
                      values=c('grey30', 'red', model.cols),
                      labels=c('M_M', 'M_M1', model.labels))+
    geom_hline(yintercept = 0, lty=2)+
    geom_vline(xintercept = 0, lty=2)+
    xlab(expression('P effect on ' * Delta * C[veg] * "  (g C " * m^-2 * ' ' * yr^-1 * ')'))+
    ylab(expression(CO[2] * '  effect on ' * Delta * C[veg] * "  (g C " * m^-2 * ' ' * yr^-1 * ')'))+
    annotate('text', x=-195, y=320, label='(h)', size=10)+
    annotate('text', x=330, y=315, label=substitute(~R^2~"= 0.67,"~italic(p)~"< 0.01"), size=9, family = 'serif')+
    theme_linedraw() +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          plot.margin = unit(c(0.8,1,0,0.5), 'cm'),
          axis.ticks.length = unit(0.2, 'cm'),
          plot.title = element_blank(),
          axis.text.x=element_text(size=22, vjust = 0),
          axis.title.x=element_text(size=22, vjust = 4),
          axis.text.y=element_text(size=22),
          axis.title.y=element_text(size=22),
          legend.text=element_text(size=22),
          legend.title=element_blank(),
          legend.position='none',
          legend.box = 'horizontal',
          legend.box.just = 'left',
          legend.key.size = unit(1.8, 'cm'))
  
  # NEP
  summary(lm(NEP_CO2~NEP_P, data = subset(P_CO2_plot, Value == 'magnitude')))
  p3 <-
    ggplot(data = subset(P_CO2_plot, Value == 'magnitude'),
           aes(NEP_P, NEP_CO2))+
    geom_point(aes(fill=ModName), pch=21, size=8)+
    geom_errorbar(aes(x=NEP_P,
                      ymax = NEP_CO2 + NEP_sd_CO2,
                      ymin = NEP_CO2 - NEP_sd_CO2),
                  linewidth = 1.2,
                  col='grey40',
                  alpha=0.5,
                  width=1.5)+
    geom_errorbar(aes(y=NEP_CO2,
                      xmax = NEP_P + NEP_sd_P,
                      xmin = NEP_P - NEP_sd_P),
                  linewidth = 1.2,
                  col='grey40',
                  alpha=0.5,
                  width=1.5)+
    geom_smooth(method = 'lm', lty=1, col='blue', linewidth=2.5, se=T, alpha=0.2)+
    scale_x_continuous(limits = c(-210,610),
                       breaks = seq(-200,600,200))+
    scale_y_continuous(limits = c(-40,340),
                       breaks = seq(0,300,100))+
    scale_fill_manual(name="",
                      values=c('grey30', 'red', model.cols),
                      labels=c('M_M', 'M_M1', model.labels))+
    geom_hline(yintercept = 0, lty=2)+
    geom_vline(xintercept = 0, lty=2)+
    xlab(expression('P effect on NEP (g C ' * m^-2 * ' ' * yr^-1 * ')'))+
    ylab(expression(CO[2] * '  effect on NEP (g C ' * m^-2 * ' ' * yr^-1 * ')'))+
    annotate('text', x=-195, y=320,  label='(i)', size=10)+
    annotate('text', x=330, y=315, label=substitute(~R^2~"= 0.80,"~italic(p)~"< 0.001"), size=9, family = 'serif')+
    theme_linedraw() +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          plot.margin = unit(c(0.8,1,0,0.5), 'cm'),
          axis.ticks.length = unit(0.2, 'cm'),
          plot.title = element_blank(),
          axis.text.x=element_text(size=22, vjust = 0),
          axis.title.x=element_text(size=22, vjust = 4),
          axis.text.y=element_text(size=22),
          axis.title.y=element_text(size=22),
          legend.text=element_text(size=22),
          legend.title=element_blank(),
          legend.position='none',
          legend.box = 'horizontal',
          legend.box.just = 'left',
          legend.key.size = unit(1.8, 'cm'))
  
  # combine plot
  plot_output1 <- 
    pc1 + pi1 + p1 +
    pc2 + pi2 + p2 + 
    pc3 + pi3 + p3 + 
    plot_layout(ncol = 3,
                widths = c(1.5,3.2,2.5,1.5,3.2,2.5,1.5,3.2,2.5))
  
}


## plot model-specific interactive effects of CO2 and HP on C and P processes
# ----------------------------------------------------------------------------------
{
  P_x_CO2_plot <- subset(PE_effect_plot, Effect == 'HIP_x_CO2' & Value == 'percentage',
                         select = c('ModName',
                                    'GPP', 'CUE', 'NPP', 'RAU', 
                                    'ABR', 'CPL', 'CPW', 'CPFR',  
                                    'RHET', 'RSOIL', 'CSOIL',
                                    paste(c('GPP', 'CUE', 'NPP', 'RAU', 
                                            'ABR', 'CPL', 'CPW', 'CPFR',  
                                            'RHET', 'RSOIL', 'CSOIL'), 'sd', sep = '_')))  

  P_x_CO2_plot1 <- melt(P_x_CO2_plot[1:12], by='ModName') 
  names(P_x_CO2_plot1)[3] <- 'mean'
  
  P_x_CO2_plot2 <- P_x_CO2_plot[c(1, 13:23)]
  names(P_x_CO2_plot2)[2:12] <- names(P_x_CO2_plot)[2:12]
  P_x_CO2_plot2 <- melt(P_x_CO2_plot2, by='ModName') 
  names(P_x_CO2_plot2)[3] <- 'sd'
  
  P_x_CO2_plot <- merge(P_x_CO2_plot1, P_x_CO2_plot2, by = c('ModName', 'variable'))
  
  # add sign
  P_x_CO2_plot$sign[P_x_CO2_plot$mean < 0] <- 'negative'
  P_x_CO2_plot$sign[P_x_CO2_plot$mean >= 0] <- 'positive'
  
  P_x_CO2_plot$variable <- factor(P_x_CO2_plot$variable,
                                  levels = rev(c('GPP', 'CUE', 'NPP', 'RAU', 
                                                 'ABR', 'CPL', 'CPW', 'CPFR',  
                                                 'RHET', 'RSOIL', 'CSOIL')),
                                  labels = rev(c('GPP', 'CUE', 'NPP', 'RAU', 
                                                 'ABR', 'CPL', 'CPW', 'CPFR',  
                                                 'RHET', 'RSOIL', 'CSOIL')))
  
  plot_m1 <-
    ggplot(data = subset(P_x_CO2_plot, ModName == 'A_GDAYP'))+
    geom_bar(stat = "identity",
             mapping=aes(variable, mean, fill=sign),
             position=position_dodge(preserve = 'single'),
             col="black",
             width = 0.5,
             alpha=0.7)+
    geom_errorbar(aes(variable,
                      ymax = mean + sd,
                      ymin = mean - sd),
                  width = 0.3, linewidth=1, col='black')+
    scale_fill_manual(names(''),
                      values = c('blue', 'red'))+
    scale_y_continuous(limits = c(-33,33),
                       breaks = seq(-30,30,15))+
    geom_hline(yintercept = 0, linetype=1, col='black')+
    geom_vline(xintercept = c(5.5, 8.5), linetype=2, linewidth=1, col='grey70')+
    annotate('text', x=10.8, y=-30, label='(a)', size=10)+
    ylab(NULL)+
    xlab(NULL)+
    ggtitle('GDAYP')+
    coord_flip()+
    theme_linedraw() +
    theme(panel.grid.major.x=element_blank(),
          panel.grid.major.y=element_blank(),
          panel.grid.minor=element_blank(),
          panel.spacing.x = unit(0.4, 'cm'),
          panel.spacing.y = unit(0.2, 'cm'),
          plot.title = element_text(size=26, vjust = 2, hjust = 0.5),
          plot.margin = unit(c(1.0,0.2,1.0,0.3), 'cm'),
          axis.text.x=element_text(size=22, vjust = -0.5),
          axis.title.x=element_text(size=26, vjust = -3),
          axis.text.y=element_text(size=21),
          # axis.text.y=element_blank(),
          # axis.text.y=element_blank(),
          axis.title.y=element_text(size=22, vjust=5),
          axis.ticks.length = unit(0.2, 'cm'),
          strip.text = element_blank(),
          legend.position = 'none')
  
  plot_m2 <-
    ggplot(data = subset(P_x_CO2_plot, ModName == 'H_QUJSM'))+
    geom_bar(stat = "identity",
             mapping=aes(variable, mean, fill=sign),
             position=position_dodge(preserve = 'single'),
             col="black",
             width = 0.5,
             alpha=0.7)+
    geom_errorbar(aes(variable,
                      ymax = mean + sd,
                      ymin = mean - sd),
                  width = 0.3, linewidth=1, col='black')+
    scale_fill_manual(names(''),
                      values = c('blue', 'red'))+
    scale_y_continuous(limits = c(-33,33),
                       breaks = seq(-30,30,15))+
    geom_hline(yintercept = 0, linetype=1, col='black')+
    geom_vline(xintercept = c(5.5, 8.5), linetype=2, linewidth=1, col='grey70')+
    annotate('text', x=10.8, y=-30, label='(b)', size=10)+
    xlab(NULL)+
    ylab(NULL)+
    ggtitle('QUJSM')+
    coord_flip()+
    theme_linedraw() +
    theme(panel.grid.major.x=element_blank(),
          panel.grid.major.y=element_blank(),
          panel.grid.minor=element_blank(),
          panel.spacing.x = unit(0.4, 'cm'),
          panel.spacing.y = unit(0.2, 'cm'),
          plot.title = element_text(size=26, vjust = 2, hjust = 0.5),
          plot.margin = unit(c(1.0,0.2,1.0,0.3), 'cm'),
          axis.text.x=element_text(size=22, vjust = -0.5),
          axis.title.x=element_text(size=26, vjust = -1.5),
          # axis.title.x=element_blank(),
          # axis.text.y=element_text(size=21),
          axis.text.y=element_blank(),
          axis.title.y=element_text(size=22, vjust=5),
          axis.ticks.length = unit(0.2, 'cm'),
          strip.text = element_blank(),
          legend.position = 'none')
  
  plot_m3 <-
    ggplot(data = subset(P_x_CO2_plot, ModName == 'C_CABLP'))+
    geom_bar(stat = "identity",
             mapping=aes(variable, mean, fill=sign),
             position=position_dodge(preserve = 'single'),
             col="black",
             width = 0.5,
             alpha=0.7)+
    geom_errorbar(aes(variable,
                      ymax = mean + sd,
                      ymin = mean - sd),
                  width = 0.3, linewidth=1, col='black')+
    scale_fill_manual(names(''),
                      values = c('blue', 'red'))+
    scale_y_continuous(limits = c(-33,33),
                       breaks = seq(-30,30,15))+
    geom_hline(yintercept = 0, linetype=1, col='black')+
    geom_vline(xintercept = c(5.5, 8.5), linetype=2, linewidth=1, col='grey70')+
    annotate('text', x=10.8, y=-30, label='(c)', size=10)+
    xlab(NULL)+
    ylab(NULL)+
    ggtitle('CABLP')+
    coord_flip()+
    theme_linedraw() +
    theme(panel.grid.major.x=element_blank(),
          panel.grid.major.y=element_blank(),
          panel.grid.minor=element_blank(),
          panel.spacing.x = unit(0.4, 'cm'),
          panel.spacing.y = unit(0.2, 'cm'),
          plot.title = element_text(size=26, vjust = 2, hjust = 0.5),
          plot.margin = unit(c(1.0,0.2,0.3,0.3), 'cm'),
          axis.text.x=element_text(size=22, vjust = -0.5),
          # axis.title.x=element_text(size=22, vjust = -1.5),
          axis.title.x=element_blank(),
          # axis.text.y=element_text(size=21),
          axis.text.y=element_blank(),
          axis.title.y=element_text(size=22, vjust=5),
          axis.ticks.length = unit(0.2, 'cm'),
          strip.text = element_blank(),
          legend.position = 'none')
  
  plot_m4 <-
    ggplot(data = subset(P_x_CO2_plot, ModName == 'E_OCHDP'))+
    geom_bar(stat = "identity",
             mapping=aes(variable, mean, fill=sign),
             position=position_dodge(preserve = 'single'),
             col="black",
             width = 0.5,
             alpha=0.7)+
    geom_errorbar(aes(variable,
                      ymax = mean + sd,
                      ymin = mean - sd),
                  width = 0.3, linewidth=1, col='black')+
    scale_fill_manual(names(''),
                      values = c('blue', 'red'))+
    scale_y_continuous(limits = c(-42,38),
                       breaks = seq(-30,30,15))+
    geom_hline(yintercept = 0, linetype=1, col='black')+
    geom_vline(xintercept = c(5.5, 8.5), linetype=2, linewidth=1, col='grey70')+
    annotate('text', x=10.8, y=-38, label='(d)', size=10)+
    xlab(NULL)+
    ylab(NULL)+
    ggtitle('OCHDP')+
    coord_flip()+
    theme_linedraw() +
    theme(panel.grid.major.x=element_blank(),
          panel.grid.major.y=element_blank(),
          panel.grid.minor=element_blank(),
          panel.spacing.x = unit(0.4, 'cm'),
          panel.spacing.y = unit(0.2, 'cm'),
          plot.title = element_text(size=26, vjust = 2, hjust = 0.5),
          plot.margin = unit(c(1.0,0.2,0.3,0.3), 'cm'),
          axis.text.x=element_text(size=22, vjust = -0.5),
          # axis.title.x=element_text(size=22, vjust = -1.5),
          axis.title.x=element_blank(),
          # axis.text.y=element_text(size=21),
          axis.text.y=element_blank(),
          # axis.text.y=element_blank(),
          axis.title.y=element_text(size=22, vjust=5),
          axis.ticks.length = unit(0.2, 'cm'),
          strip.text = element_blank(),
          legend.position = 'none')
  
  plot_m5 <-
    ggplot(data = subset(P_x_CO2_plot, ModName == 'G_OCHDX'))+
    geom_bar(stat = "identity",
             mapping=aes(variable, mean, fill=sign),
             position=position_dodge(preserve = 'single'),
             col="black",
             width = 0.5,
             alpha=0.7)+
    geom_errorbar(aes(variable,
                      ymax = mean + sd,
                      ymin = mean - sd),
                  width = 0.3, linewidth=1, col='black')+
    scale_fill_manual(names(''),
                      values = c('blue', 'red'))+
    scale_y_continuous(limits = c(-33,33),
                       breaks = seq(-30,30,15))+
    geom_hline(yintercept = 0, linetype=1, col='black')+
    geom_vline(xintercept = c(5.5, 8.5), linetype=2, linewidth=1, col='grey70')+
    annotate('text', x=10.8, y=-30, label='(e)', size=10)+
    xlab(NULL)+
    ylab('Effect size (%)')+
    ggtitle('OCHDX')+
    coord_flip()+
    theme_linedraw() +
    theme(panel.grid.major.x=element_blank(),
          panel.grid.major.y=element_blank(),
          panel.grid.minor=element_blank(),
          panel.spacing.x = unit(0.4, 'cm'),
          panel.spacing.y = unit(0.2, 'cm'),
          plot.title = element_text(size=26, vjust = 2, hjust = 0.5),
          plot.margin = unit(c(1.0,0.2,1.0,0.3), 'cm'),
          axis.text.x=element_text(size=22, vjust = -0.5),
          axis.title.x=element_text(size=26, vjust = -3),
          axis.text.y=element_text(size=21),
          # axis.text.y=element_blank(),
          # axis.text.y=element_blank(),
          axis.title.y=element_text(size=22, vjust=5),
          axis.ticks.length = unit(0.2, 'cm'),
          strip.text = element_blank(),
          legend.position = 'none')
  
  plot_m6 <-
    ggplot(data = subset(P_x_CO2_plot, ModName == 'B_ELMV1'))+
    geom_bar(stat = "identity",
             mapping=aes(variable, mean, fill=sign),
             position=position_dodge(preserve = 'single'),
             col="black",
             width = 0.5,
             alpha=0.7)+
    geom_errorbar(aes(variable,
                      ymax = mean + sd,
                      ymin = mean - sd),
                  width = 0.3, linewidth=1, col='black')+
    scale_fill_manual(names(''),
                      values = c('blue', 'red'))+
    scale_y_continuous(limits = c(-33,33),
                       breaks = seq(-30,30,15))+
    geom_hline(yintercept = 0, linetype=1, col='black')+
    geom_vline(xintercept = c(5.5, 8.5), linetype=2, linewidth=1, col='grey70')+
    annotate('text', x=10.8, y=-30, label='(f)', size=10)+
    ylab('Effect size (%)')+
    xlab(NULL)+
    ggtitle('ELMV1')+
    coord_flip()+
    theme_linedraw() +
    theme(panel.grid.major.x=element_blank(),
          panel.grid.major.y=element_blank(),
          panel.grid.minor=element_blank(),
          panel.spacing.x = unit(0.4, 'cm'),
          panel.spacing.y = unit(0.2, 'cm'),
          plot.title = element_text(size=26, vjust = 2, hjust = 0.5),
          plot.margin = unit(c(1.0,0.2,1.0,0.3), 'cm'),
          axis.text.x=element_text(size=22, vjust = -0.5),
          axis.title.x=element_text(size=26, vjust = -3),
          # axis.text.y=element_text(size=21),
          axis.text.y=element_blank(),
          # axis.text.y=element_blank(),
          axis.title.y=element_text(size=22, vjust=5),
          axis.ticks.length = unit(0.2, 'cm'),
          strip.text = element_blank(),
          legend.position = 'none')
  
  plot_m7 <-
    ggplot(data = subset(P_x_CO2_plot, ModName == 'F_QUINC'))+
    geom_bar(stat = "identity",
             mapping=aes(variable, mean, fill=sign),
             position=position_dodge(preserve = 'single'),
             col="black",
             width = 0.5,
             alpha=0.7)+
    geom_errorbar(aes(variable,
                      ymax = mean + sd,
                      ymin = mean - sd),
                  width = 0.3, linewidth=1, col='black')+
    scale_fill_manual(names(''),
                      values = c('blue', 'red'))+
    scale_y_continuous(limits = c(-53,28),
                       breaks = seq(-50,25,25))+
    geom_hline(yintercept = 0, linetype=1, col='black')+
    geom_vline(xintercept = c(5.5, 8.5), linetype=2, linewidth=1, col='grey70')+
    annotate('text', x=10.8, y=-49, label='(g)', size=10)+
    ylab('Effect size (%)')+
    xlab(NULL)+
    ggtitle('QUINC')+
    coord_flip()+
    theme_linedraw() +
    theme(panel.grid.major.x=element_blank(),
          panel.grid.major.y=element_blank(),
          panel.grid.minor=element_blank(),
          panel.spacing.x = unit(0.4, 'cm'),
          panel.spacing.y = unit(0.2, 'cm'),
          plot.title = element_text(size=26, vjust = 2, hjust = 0.5),
          plot.margin = unit(c(1.0,1.0,1.0,0.3), 'cm'),
          axis.text.x=element_text(size=22, vjust = -0.5),
          axis.title.x=element_text(size=26, vjust = -3),
          # axis.text.y=element_text(size=21),
          axis.text.y=element_blank(),
          # axis.text.y=element_blank(),
          axis.title.y=element_text(size=22, vjust=5),
          axis.ticks.length = unit(0.2, 'cm'),
          strip.text = element_blank(),
          legend.position = 'none')
  
  plot_m8 <-
    ggplot(data = subset(P_x_CO2_plot, ModName == 'D_LPJGP'))+
    geom_bar(stat = "identity",
             mapping=aes(variable, mean, fill=sign),
             position=position_dodge(preserve = 'single'),
             col="black",
             width = 0.5,
             alpha=0.7)+
    geom_errorbar(aes(variable,
                      ymax = mean + sd,
                      ymin = mean - sd),
                  width = 0.3, linewidth=1, col='black')+
    scale_fill_manual(names(''),
                      values = c('blue', 'red'))+
    scale_y_continuous(limits = c(-305,105),
                       breaks = c(-300,-200,-100, 0, 100))+
    geom_hline(yintercept = 0, linetype=1, col='black')+
    geom_vline(xintercept = c(5.5, 8.5), linetype=2, linewidth=1, col='grey70')+
    annotate('text', x=10.8, y=-300, label='(h)', size=10)+
    ylab('Effect size (%)')+
    xlab(NULL)+
    ggtitle('LPJGP')+
    coord_flip()+
    theme_linedraw() +
    theme(panel.grid.major.x=element_blank(),
          panel.grid.major.y=element_blank(),
          panel.grid.minor=element_blank(),
          panel.spacing.x = unit(0.4, 'cm'),
          panel.spacing.y = unit(0.2, 'cm'),
          plot.title = element_text(size=26, vjust = 2, hjust = 0.5),
          plot.margin = unit(c(1.0,0.2,1.0,0.3), 'cm'),
          axis.text.x=element_text(size=22, vjust = -0.5),
          axis.title.x=element_text(size=26, vjust = -3),
          # axis.text.y=element_text(size=21),
          axis.text.y=element_blank(),
          axis.title.y=element_text(size=22, vjust=5),
          axis.ticks.length = unit(0.2, 'cm'),
          strip.text = element_blank(),
          legend.position = 'none')
  
  # combine plot
  plot_interaction <- 
    plot_m1 + plot_m2 + plot_m3 + plot_m4 +
    plot_m5 + plot_m6 + plot_m7 + plot_m8 +
    plot_layout(ncol = 4) 
}


## plot model-specific CO2 effect under ambient and enriched P
# ----------------------------------------------------------------------------------
{
  P_x_CO2_plot <- subset(PE_effect_plot, (Effect == 'CO2' | Effect == 'CO2_HIP') & Value == 'percentage',
                         select = c('ModName', 'Effect',
                                    'GPP', 'NPP', 'RAU', 'CUE', 'RHET',
                                    'LAI', 'Cveg', 'ABR', 'CSOIL',
                                    'CPL', 'CNL', 'CPFR',  
                                    'PUP', 'NUP', 'PLAB', 'PMIN',
                                    paste(c('GPP', 'NPP', 'RAU', 'CUE', 'RHET',
                                            'LAI', 'Cveg', 'ABR', 'CSOIL',
                                            'CPL', 'CNL', 'CPFR',  
                                            'PUP', 'NUP', 'PLAB', 'PMIN'), 'sd', sep = '_')))  
  
  P_x_CO2_plot1 <- melt(P_x_CO2_plot[1:18], by=c('ModName', 'Effect')) 
  names(P_x_CO2_plot1)[4] <- 'mean'
  
  P_x_CO2_plot2 <- P_x_CO2_plot[c(1:2, 19:34)]
  names(P_x_CO2_plot2)[3:18] <- names(P_x_CO2_plot)[3:18]
  P_x_CO2_plot2 <- melt(P_x_CO2_plot2, by=c('ModName', 'Effect')) 
  names(P_x_CO2_plot2)[4] <- 'sd'
  
  P_x_CO2_plot <- merge(P_x_CO2_plot1, P_x_CO2_plot2, by = c('ModName', 'Effect', 'variable'))
  
  # add sign
  P_x_CO2_plot$sign[P_x_CO2_plot$mean < 0] <- 'negative'
  P_x_CO2_plot$sign[P_x_CO2_plot$mean >= 0] <- 'positive'
  
  P_x_CO2_plot$variable <- factor(P_x_CO2_plot$variable,
                                  levels = rev(c('GPP', 'NPP', 'RAU', 'CUE', 'RHET',
                                                 'LAI', 'Cveg', 'ABR', 'CSOIL',
                                                 'CPL', 'CNL', 'CPFR',  
                                                 'PUP', 'NUP', 'PLAB', 'PMIN')),
                                  labels = rev(c('GPP', 'NPP', 'RAU', 'CUE', 'RHET',
                                                 'LAI', 'Cveg', 'ABR', 'CSOIL',
                                                 'CPL', 'CNL', 'CPFR',  
                                                 'PUP', 'NUP', 'PLAB', 'PMIN')))
  P_x_CO2_plot$Type <- factor(paste(P_x_CO2_plot$Effect, P_x_CO2_plot$sign, sep = '_'),
                              levels = c('CO2_negative', 'CO2_positive',
                                         'CO2_HIP_negative', 'CO2_HIP_positive'))
  
  plot_m1 <-
    ggplot(data = subset(P_x_CO2_plot, ModName == 'A_GDAYP'))+
    geom_bar(stat = "identity",
             mapping=aes(variable, mean, group=Type, fill=Type),
             position=position_dodge(preserve = 'single'),
             col="black",
             width = 0.5,
             alpha=0.7)+
    geom_errorbar(aes(variable,
                      ymax = mean + sd,
                      ymin = mean - sd,
                      group = Type),
                  position=position_dodge(preserve = 'single', width = 0.5),
                  width = 0.3, col='black')+
    scale_fill_manual(names(''),
                      values = c('white', 'white', 'blue', 'red'))+
    scale_y_continuous(limits = c(-33,33),
                       breaks = seq(-30,30,15))+
    geom_hline(yintercept = 0, linetype=1, col='black')+
    geom_vline(xintercept = c(5.5, 8.5, 12.5), linetype=2, linewidth=1, col='grey70')+
    annotate('text', x=15.5, y=-30, label='(a)', size=10)+
    ylab(NULL)+
    xlab(NULL)+
    ggtitle('GDAYP')+
    coord_flip()+
    theme_linedraw() +
    theme(panel.grid.major.x=element_blank(),
          panel.grid.major.y=element_blank(),
          panel.grid.minor=element_blank(),
          panel.spacing.x = unit(0.4, 'cm'),
          panel.spacing.y = unit(0.2, 'cm'),
          plot.title = element_text(size=26, vjust = 2, hjust = 0.5),
          plot.margin = unit(c(1.0,0.2,1.0,0.3), 'cm'),
          axis.text.x=element_text(size=22, vjust = -0.5),
          axis.title.x=element_text(size=26, vjust = -3),
          axis.text.y=element_text(size=21),
          # axis.text.y=element_blank(),
          # axis.text.y=element_blank(),
          axis.title.y=element_text(size=22, vjust=5),
          axis.ticks.length = unit(0.2, 'cm'),
          strip.text = element_blank(),
          legend.position = 'none')
  
  plot_m2 <-
    ggplot(data = subset(P_x_CO2_plot, ModName == 'H_QUJSM'))+
    geom_bar(stat = "identity",
             mapping=aes(variable, mean, group=Type, fill=Type),
             position=position_dodge(preserve = 'single'),
             col="black",
             width = 0.5,
             alpha=0.7)+
    geom_errorbar(aes(variable,
                      ymax = mean + sd,
                      ymin = mean - sd,
                      group = Type),
                  position=position_dodge(preserve = 'single', width = 0.5),
                  width = 0.3, col='black')+
    scale_fill_manual(names(''),
                      values = c('white', 'white', 'blue', 'red'))+
    scale_y_continuous(limits = c(-18,48),
                       breaks = seq(-15,45,15))+
    geom_hline(yintercept = 0, linetype=1, col='black')+
    geom_vline(xintercept = c(5.5, 8.5, 12.5), linetype=2, linewidth=1, col='grey70')+
    annotate('text', x=15.5, y=-14, label='(b)', size=10)+
    xlab(NULL)+
    ylab(NULL)+
    ggtitle('QUJSM')+
    coord_flip()+
    theme_linedraw() +
    theme(panel.grid.major.x=element_blank(),
          panel.grid.major.y=element_blank(),
          panel.grid.minor=element_blank(),
          panel.spacing.x = unit(0.4, 'cm'),
          panel.spacing.y = unit(0.2, 'cm'),
          plot.title = element_text(size=26, vjust = 2, hjust = 0.5),
          plot.margin = unit(c(1.0,0.2,1.0,0.3), 'cm'),
          axis.text.x=element_text(size=22, vjust = -0.5),
          axis.title.x=element_text(size=26, vjust = -1.5),
          # axis.title.x=element_blank(),
          # axis.text.y=element_text(size=21),
          axis.text.y=element_blank(),
          axis.title.y=element_text(size=22, vjust=5),
          axis.ticks.length = unit(0.2, 'cm'),
          strip.text = element_blank(),
          legend.position = 'none')
  
  plot_m3 <-
    ggplot(data = subset(P_x_CO2_plot, ModName == 'C_CABLP'))+
    geom_bar(stat = "identity",
             mapping=aes(variable, mean, group=Type, fill=Type),
             position=position_dodge(preserve = 'single'),
             col="black",
             width = 0.5,
             alpha=0.7)+
    geom_errorbar(aes(variable,
                      ymax = mean + sd,
                      ymin = mean - sd,
                      group = Type),
                  position=position_dodge(preserve = 'single', width = 0.5),
                  width = 0.3, col='black')+
    scale_fill_manual(names(''),
                      values = c('white', 'white', 'blue', 'red'))+
    scale_y_continuous(limits = c(-125,45),
                       breaks = seq(-120,40,40))+
    geom_hline(yintercept = 0, linetype=1, col='black')+
    geom_vline(xintercept = c(5.5, 8.5, 12.5), linetype=2, linewidth=1, col='grey70')+
    annotate('text', x=15.5, y=-115, label='(c)', size=10)+
    xlab(NULL)+
    ylab(NULL)+
    ggtitle('CABLP')+
    coord_flip()+
    theme_linedraw() +
    theme(panel.grid.major.x=element_blank(),
          panel.grid.major.y=element_blank(),
          panel.grid.minor=element_blank(),
          panel.spacing.x = unit(0.4, 'cm'),
          panel.spacing.y = unit(0.2, 'cm'),
          plot.title = element_text(size=26, vjust = 2, hjust = 0.5),
          plot.margin = unit(c(1.0,0.2,0.3,0.3), 'cm'),
          axis.text.x=element_text(size=22, vjust = -0.5),
          # axis.title.x=element_text(size=22, vjust = -1.5),
          axis.title.x=element_blank(),
          # axis.text.y=element_text(size=21),
          axis.text.y=element_blank(),
          axis.title.y=element_text(size=22, vjust=5),
          axis.ticks.length = unit(0.2, 'cm'),
          strip.text = element_blank(),
          legend.position = 'none')
  
  plot_m4 <-
    ggplot(data = subset(P_x_CO2_plot, ModName == 'E_OCHDP'))+
    geom_bar(stat = "identity",
             mapping=aes(variable, mean, group=Type, fill=Type),
             position=position_dodge(preserve = 'single'),
             col="black",
             width = 0.5,
             alpha=0.7)+
    geom_errorbar(aes(variable,
                      ymax = mean + sd,
                      ymin = mean - sd,
                      group = Type),
                  position=position_dodge(preserve = 'single', width = 0.5),
                  width = 0.3, col='black')+
    scale_fill_manual(names(''),
                      values = c('white', 'white', 'blue', 'red'))+
    scale_y_continuous(limits = c(-45,65),
                       breaks = seq(-40,60,20))+
    geom_hline(yintercept = 0, linetype=1, col='black')+
    geom_vline(xintercept = c(5.5, 8.5, 12.5), linetype=2, linewidth=1, col='grey70')+
    annotate('text', x=15.5, y=-40, label='(d)', size=10)+
    xlab(NULL)+
    ylab(NULL)+
    ggtitle('OCHDP')+
    coord_flip()+
    theme_linedraw() +
    theme(panel.grid.major.x=element_blank(),
          panel.grid.major.y=element_blank(),
          panel.grid.minor=element_blank(),
          panel.spacing.x = unit(0.4, 'cm'),
          panel.spacing.y = unit(0.2, 'cm'),
          plot.title = element_text(size=26, vjust = 2, hjust = 0.5),
          plot.margin = unit(c(1.0,0.2,0.3,0.3), 'cm'),
          axis.text.x=element_text(size=22, vjust = -0.5),
          # axis.title.x=element_text(size=22, vjust = -1.5),
          axis.title.x=element_blank(),
          # axis.text.y=element_text(size=21),
          axis.text.y=element_blank(),
          # axis.text.y=element_blank(),
          axis.title.y=element_text(size=22, vjust=5),
          axis.ticks.length = unit(0.2, 'cm'),
          strip.text = element_blank(),
          legend.position = 'none')
  
  plot_m5 <-
    ggplot(data = subset(P_x_CO2_plot, ModName == 'G_OCHDX'))+
    geom_bar(stat = "identity",
             mapping=aes(variable, mean, group=Type, fill=Type),
             position=position_dodge(preserve = 'single'),
             col="black",
             width = 0.5,
             alpha=0.7)+
    geom_errorbar(aes(variable,
                      ymax = mean + sd,
                      ymin = mean - sd,
                      group = Type),
                  position=position_dodge(preserve = 'single', width = 0.5),
                  width = 0.3, col='black')+
    scale_fill_manual(names(''),
                      values = c('white', 'white', 'blue', 'red'))+
    scale_y_continuous(limits = c(-43,43),
                       breaks = seq(-40,40,20))+
    geom_hline(yintercept = 0, linetype=1, col='black')+
    geom_vline(xintercept = c(5.5, 8.5, 12.5), linetype=2, linewidth=1, col='grey70')+
    annotate('text', x=15.5, y=-38, label='(e)', size=10)+
    xlab(NULL)+
    ylab('Effect size (%)')+
    ggtitle('OCHDX')+
    coord_flip()+
    theme_linedraw() +
    theme(panel.grid.major.x=element_blank(),
          panel.grid.major.y=element_blank(),
          panel.grid.minor=element_blank(),
          panel.spacing.x = unit(0.4, 'cm'),
          panel.spacing.y = unit(0.2, 'cm'),
          plot.title = element_text(size=26, vjust = 2, hjust = 0.5),
          plot.margin = unit(c(1.0,0.2,1.0,0.3), 'cm'),
          axis.text.x=element_text(size=22, vjust = -0.5),
          axis.title.x=element_text(size=26, vjust = -3),
          axis.text.y=element_text(size=21),
          # axis.text.y=element_blank(),
          # axis.text.y=element_blank(),
          axis.title.y=element_text(size=22, vjust=5),
          axis.ticks.length = unit(0.2, 'cm'),
          strip.text = element_blank(),
          legend.position = 'none')
  
  plot_m6 <-
    ggplot(data = subset(P_x_CO2_plot, ModName == 'B_ELMV1'))+
    geom_bar(stat = "identity",
             mapping=aes(variable, mean, group=Type, fill=Type),
             position=position_dodge(preserve = 'single'),
             col="black",
             width = 0.5,
             alpha=0.7)+
    geom_errorbar(aes(variable,
                      ymax = mean + sd,
                      ymin = mean - sd,
                      group = Type),
                  position=position_dodge(preserve = 'single', width = 0.5),
                  width = 0.3, col='black')+
    scale_fill_manual(names(''),
                      values = c('white', 'white', 'blue', 'red'))+
    scale_y_continuous(limits = c(-25,81),
                       breaks = seq(-20,80,20))+
    geom_hline(yintercept = 0, linetype=1, col='black')+
    geom_vline(xintercept = c(5.5, 8.5, 12.5), linetype=2, linewidth=1, col='grey70')+
    annotate('text', x=15.5, y=-20, label='(f)', size=10)+
    ylab('Effect size (%)')+
    xlab(NULL)+
    ggtitle('ELMV1')+
    coord_flip()+
    theme_linedraw() +
    theme(panel.grid.major.x=element_blank(),
          panel.grid.major.y=element_blank(),
          panel.grid.minor=element_blank(),
          panel.spacing.x = unit(0.4, 'cm'),
          panel.spacing.y = unit(0.2, 'cm'),
          plot.title = element_text(size=26, vjust = 2, hjust = 0.5),
          plot.margin = unit(c(1.0,0.2,1.0,0.3), 'cm'),
          axis.text.x=element_text(size=22, vjust = -0.5),
          axis.title.x=element_text(size=26, vjust = -3),
          # axis.text.y=element_text(size=21),
          axis.text.y=element_blank(),
          # axis.text.y=element_blank(),
          axis.title.y=element_text(size=22, vjust=5),
          axis.ticks.length = unit(0.2, 'cm'),
          strip.text = element_blank(),
          legend.position = 'none')
  
  plot_m7 <-
    ggplot(data = subset(P_x_CO2_plot, ModName == 'F_QUINC'))+
    geom_bar(stat = "identity",
             mapping=aes(variable, mean, group=Type, fill=Type),
             position=position_dodge(preserve = 'single'),
             col="black",
             width = 0.5,
             alpha=0.7)+
    geom_errorbar(aes(variable,
                      ymax = mean + sd,
                      ymin = mean - sd,
                      group = Type),
                  position=position_dodge(preserve = 'single', width = 0.5),
                  width = 0.3, col='black')+
    scale_fill_manual(names(''),
                      values = c('white', 'white', 'blue', 'red'))+
    scale_y_continuous(limits = c(-43,43),
                       breaks = seq(-40,40,20))+
    geom_hline(yintercept = 0, linetype=1, col='black')+
    geom_vline(xintercept = c(5.5, 8.5, 12.5), linetype=2, linewidth=1, col='grey70')+
    annotate('text', x=15.5, y=-38, label='(g)', size=10)+
    ylab('Effect size (%)')+
    xlab(NULL)+
    ggtitle('QUINC')+
    coord_flip()+
    theme_linedraw() +
    theme(panel.grid.major.x=element_blank(),
          panel.grid.major.y=element_blank(),
          panel.grid.minor=element_blank(),
          panel.spacing.x = unit(0.4, 'cm'),
          panel.spacing.y = unit(0.2, 'cm'),
          plot.title = element_text(size=26, vjust = 2, hjust = 0.5),
          plot.margin = unit(c(1.0,1.0,1.0,0.3), 'cm'),
          axis.text.x=element_text(size=22, vjust = -0.5),
          axis.title.x=element_text(size=26, vjust = -3),
          # axis.text.y=element_text(size=21),
          axis.text.y=element_blank(),
          # axis.text.y=element_blank(),
          axis.title.y=element_text(size=22, vjust=5),
          axis.ticks.length = unit(0.2, 'cm'),
          strip.text = element_blank(),
          legend.position = 'none')
  
  plot_m8 <-
    ggplot(data = subset(P_x_CO2_plot, ModName == 'D_LPJGP'))+
    geom_bar(stat = "identity",
             mapping=aes(variable, mean, group=Type, fill=Type),
             position=position_dodge(preserve = 'single'),
             col="black",
             width = 0.5,
             alpha=0.7)+
    geom_errorbar(aes(variable,
                      ymax = mean + sd,
                      ymin = mean - sd,
                      group = Type),
                  position=position_dodge(preserve = 'single', width = 0.5),
                  width = 0.3, col='black')+
    scale_fill_manual(names(''),
                      values = c('white', 'white', 'blue', 'red'))+
    scale_y_continuous(limits = c(-105,305),
                       breaks = c(-100, 0, 100, 200, 300))+
    geom_hline(yintercept = 0, linetype=1, col='black')+
    geom_vline(xintercept = c(5.5, 8.5, 12.5), linetype=2, linewidth=1, col='grey70')+
    annotate('text', x=15.5, y=-92, label='(h)', size=10)+
    ylab('Effect size (%)')+
    xlab(NULL)+
    ggtitle('LPJGP')+
    coord_flip()+
    theme_linedraw() +
    theme(panel.grid.major.x=element_blank(),
          panel.grid.major.y=element_blank(),
          panel.grid.minor=element_blank(),
          panel.spacing.x = unit(0.4, 'cm'),
          panel.spacing.y = unit(0.2, 'cm'),
          plot.title = element_text(size=26, vjust = 2, hjust = 0.5),
          plot.margin = unit(c(1.0,0.2,1.0,0.3), 'cm'),
          axis.text.x=element_text(size=22, vjust = -0.5),
          axis.title.x=element_text(size=26, vjust = -3),
          # axis.text.y=element_text(size=21),
          axis.text.y=element_blank(),
          axis.title.y=element_text(size=22, vjust=5),
          axis.ticks.length = unit(0.2, 'cm'),
          strip.text = element_blank(),
          legend.position = 'none')
  
  # combine plot
  plot_interaction <- 
    plot_m1 + plot_m2 + plot_m3 + plot_m4 +
    plot_m5 + plot_m6 + plot_m7 + plot_m8 +
    plot_layout(ncol = 4)  
}
  

## plot interactive effects of CO2 and P (MP/HP and 2020-2022/2023-2029) 
# ----------------------------------------------------------------------------------
{
 PE_effect_plot <- PE_effect_all

 # MDP, 2020-2022
 # ---------------------------------------------
 # GPP
 pi1 <-
  ggplot(data = subset(PE_effect_plot, Period == 'temp' & Effect == 'MDP_x_CO2' & Value == 'percentage'))+
  geom_bar(stat = "identity",
           mapping=aes(ModName, GPP, fill=ModName),
           position=position_dodge(preserve = 'single'),
           col="black",
           width = 0.6) +
  geom_errorbar(mapping=aes(x = ModName, 
                            ymax = GPP + GPP_sd,
                            ymin = GPP - GPP_sd,
                            col= 'black'),
                position=position_dodge(preserve = 'single'),
                col="black",
                width = 0.4)+
  geom_vline(xintercept=c(2.5), lty=3, linewidth=1)+
  scale_fill_manual(name="",
                    values=c('grey30', 'red', model.cols),
                    labels=c('M_M', 'M_M1', model.labels))+
  scale_x_discrete(name="",
                   # guide = guide_axis(n.dodge = 2),
                   labels=c("A_GDAYP" = "GDAYP",
                            "B_ELMV1" = "ELMV1",
                            "C_CABLP" = "CABLP",
                            "D_LPJGP" = "LPJGP",
                            "E_OCHDP" = "OCHDP",
                            "F_QUINC" = "QUINC",
                            "G_OCHDX" = "OCHDX",
                            "H_QUJSM" = "QUJSM",
                            "M_M" = "M_M",
                            "M_M1" = 'W_M'))+
  scale_y_continuous(limits = c(-20, 10),
                     breaks = seq(-20,10,10))+
  geom_hline(yintercept=0, lty=1)+
  ylab(NULL)+
  ylab(expression("Interactive effect on " * GPP * " (%)"))+
  ggtitle('MP, 2020-2022') +
  annotate('text', x=1, y=9, label='(a)', size=10)+
  theme_linedraw() +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.margin = unit(c(1,0.2,0,1.5), 'cm'),
        plot.title = element_text(size=30,
                                  hjust = 0.5),
        # axis.text.x=element_text(size=20, angle = 35, vjust = 1.0, hjust = 1.0),
        axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=22),
        axis.title.y=element_text(size=22, vjust = 1),
        axis.ticks.length = unit(0.2, 'cm'),
        legend.text=element_text(size=23),
        legend.title=element_text(size=20),
        legend.position='none',
        legend.box = 'horizontal',
        legend.box.just = 'left',
        legend.key.size = unit(1.3, 'cm'))

 # biomass growth
 pi2 <-
  ggplot(data = subset(PE_effect_plot, Period == 'temp' & Effect == 'MDP_x_CO2' & Value == 'magnitude'))+
  geom_bar(stat = "identity",
           mapping=aes(ModName, deltaCveg/1000, fill=ModName),
           position=position_dodge(preserve = 'single'),
           col="black",
           width = 0.6) +
  geom_errorbar(mapping=aes(x = ModName, 
                            ymax = deltaCveg/1000 + deltaCveg_sd/1000,
                            ymin = deltaCveg/1000 - deltaCveg_sd/1000,
                            col= 'black'),
                position=position_dodge(preserve = 'single'),
                col="black",
                width = 0.4)+
  geom_vline(xintercept=c(2.5), lty=3, linewidth=1)+
  scale_fill_manual(name="",
                    values=c('grey30', 'red', model.cols),
                    labels=c('M_M', 'M_M1', model.labels))+
  scale_x_discrete(name="",
                   # guide = guide_axis(n.dodge = 2),
                   labels=c("A_GDAYP" = "GDAYP",
                            "B_ELMV1" = "ELMV1",
                            "C_CABLP" = "CABLP",
                            "D_LPJGP" = "LPJGP",
                            "E_OCHDP" = "OCHDP",
                            "F_QUINC" = "QUINC",
                            "G_OCHDX" = "OCHDX",
                            "H_QUJSM" = "QUJSM",
                            "M_M" = "M_M",
                            "M_M1" = 'W_M'))+
  scale_y_continuous(limits = c(-0.31, 0.31),
                     breaks = seq(-0.3,0.3,0.1),
                     labels = c(-0.3,-0.2,-0.1,0,0.1,0.2,0.3))+
  geom_hline(yintercept=0, lty=1)+
  ylab(NULL)+
  ylab(expression("Interactive effect on " * Delta * C[veg] * " (kg C " * m^-2 * " " * yr^-1 * ")"))+
  annotate('text', x=1, y=0.28, label='(b)', size=10)+
  theme_linedraw() +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.margin = unit(c(0.5,0.2,0,1.5), 'cm'),
        plot.title = element_text(size=24, face="bold",
                                  hjust = 0.5),
        # axis.text.x=element_text(size=20, angle = 35, vjust = 1.0, hjust = 1.0),
        axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=22),
        axis.title.y=element_text(size=22, vjust = 1),
        axis.ticks.length = unit(0.2, 'cm'),
        legend.text=element_text(size=23),
        legend.title=element_text(size=20),
        legend.position='none',
        legend.box = 'horizontal',
        legend.box.just = 'left',
        legend.key.size = unit(1.3, 'cm'))

 # NEP
 pi3 <-
  ggplot(data = subset(PE_effect_plot, Period == 'temp' & Effect == 'MDP_x_CO2' & Value == 'magnitude'))+
  geom_bar(stat = "identity",
           mapping=aes(ModName, NEP/1000, fill=ModName),
           position=position_dodge(preserve = 'single'),
           col="black",
           width = 0.6) +
  geom_errorbar(mapping=aes(x = ModName, 
                            ymax = NEP/1000 + NEP_sd/1000,
                            ymin = NEP/1000 - NEP_sd/1000,
                            col= 'black'),
                position=position_dodge(preserve = 'single'),
                col="black",
                width = 0.4)+
  geom_vline(xintercept=c(2.5), lty=3, linewidth=1)+
  scale_fill_manual(name="",
                    values=c('grey30', 'red', model.cols),
                    labels=c('M_M', 'M_M1', model.labels))+
  scale_x_discrete(name="",
                   # guide = guide_axis(n.dodge = 2),
                   labels=c("A_GDAYP" = "GDAYP",
                            "B_ELMV1" = "ELMV1",
                            "C_CABLP" = "CABLP",
                            "D_LPJGP" = "LPJGP",
                            "E_OCHDP" = "OCHDP",
                            "F_QUINC" = "QUINC",
                            "G_OCHDX" = "OCHDX",
                            "H_QUJSM" = "QUJSM",
                            "M_M" = "M_M",
                            "M_M1" = 'W_M'))+
  scale_y_continuous(limits = c(-0.41, 0.41),
                     breaks = seq(-0.4,0.4,0.2))+
  geom_hline(yintercept=0, lty=1)+
  ylab(NULL)+
  ylab(expression("Interactive effect on " * NEP * " (kg C " * m^-2 * " " * yr^-1 * ")"))+
  annotate('text', x=1, y=0.37, label='(c)', size=10)+
  theme_linedraw() +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.margin = unit(c(0.5,0.2,0.2,1.5), 'cm'),
        plot.title = element_text(size=24, face="bold",
                                  hjust = 0.5),
        axis.text.x=element_text(size=20, angle = 35, vjust = 1.0, hjust = 1.0),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=22),
        axis.title.y=element_text(size=22, vjust = 1),
        axis.ticks.length = unit(0.2, 'cm'),
        legend.text=element_text(size=23),
        legend.title=element_text(size=20),
        legend.position='none',
        legend.box = 'horizontal',
        legend.box.just = 'left',
        legend.key.size = unit(1.3, 'cm'))


 # MDP, 2023-2029
 # -----------------------------------------------
 # GPP
 pi4 <-
  ggplot(data = subset(PE_effect_plot, Period == 'futu' & Effect == 'MDP_x_CO2' & Value == 'percentage'))+
  geom_bar(stat = "identity",
           mapping=aes(ModName, GPP, fill=ModName),
           position=position_dodge(preserve = 'single'),
           col="black",
           width = 0.6) +
  geom_errorbar(mapping=aes(x = ModName, 
                            ymax = GPP + GPP_sd,
                            ymin = GPP - GPP_sd,
                            col= 'black'),
                position=position_dodge(preserve = 'single'),
                col="black",
                width = 0.4)+
  geom_vline(xintercept=c(2.5), lty=3, linewidth=1)+
  scale_fill_manual(name="",
                    values=c('grey30', 'red', model.cols),
                    labels=c('M_M', 'M_M1', model.labels))+
  scale_x_discrete(name="",
                   # guide = guide_axis(n.dodge = 2),
                   labels=c("A_GDAYP" = "GDAYP",
                            "B_ELMV1" = "ELMV1",
                            "C_CABLP" = "CABLP",
                            "D_LPJGP" = "LPJGP",
                            "E_OCHDP" = "OCHDP",
                            "F_QUINC" = "QUINC",
                            "G_OCHDX" = "OCHDX",
                            "H_QUJSM" = "QUJSM",
                            "M_M" = "M_M",
                            "M_M1" = 'W_M'))+
  scale_y_continuous(limits = c(-20, 10),
                     breaks = seq(-20,10,10))+
  geom_hline(yintercept=0, lty=1)+
  ylab(NULL)+
  ylab(expression("Interactive effect on " * GPP * " (%)"))+
  ggtitle('MP, 2023-2029') +
  annotate('text', x=1, y=9, label='(d)', size=10)+
  theme_linedraw() +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.margin = unit(c(1,0.2,0,0.2), 'cm'),
        plot.title = element_text(size=30, 
                                  hjust = 0.5),
        # axis.text.x=element_text(size=20, angle = 35, vjust = 1.0, hjust = 1.0),
        axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=22),
        # axis.title.y=element_text(size=22, vjust = 1),
        axis.title.y = element_blank(),
        axis.ticks.length = unit(0.2, 'cm'),
        legend.text=element_text(size=23),
        legend.title=element_text(size=20),
        legend.position='none',
        legend.box = 'horizontal',
        legend.box.just = 'left',
        legend.key.size = unit(1.3, 'cm'))

 # biomass growth
 pi5 <-
  ggplot(data = subset(PE_effect_plot, Period == 'futu' & Effect == 'MDP_x_CO2' & Value == 'magnitude'))+
  geom_bar(stat = "identity",
           mapping=aes(ModName, deltaCveg/1000, fill=ModName),
           position=position_dodge(preserve = 'single'),
           col="black",
           width = 0.6) +
  geom_errorbar(mapping=aes(x = ModName, 
                            ymax = deltaCveg/1000 + deltaCveg_sd/1000,
                            ymin = deltaCveg/1000 - deltaCveg_sd/1000,
                            col= 'black'),
                position=position_dodge(preserve = 'single'),
                col="black",
                width = 0.4)+
  geom_vline(xintercept=c(2.5), lty=3, linewidth=1)+
  scale_fill_manual(name="",
                    values=c('grey30', 'red', model.cols),
                    labels=c('M_M', 'M_M1', model.labels))+
  scale_x_discrete(name="",
                   # guide = guide_axis(n.dodge = 2),
                   labels=c("A_GDAYP" = "GDAYP",
                            "B_ELMV1" = "ELMV1",
                            "C_CABLP" = "CABLP",
                            "D_LPJGP" = "LPJGP",
                            "E_OCHDP" = "OCHDP",
                            "F_QUINC" = "QUINC",
                            "G_OCHDX" = "OCHDX",
                            "H_QUJSM" = "QUJSM",
                            "M_M" = "M_M",
                            "M_M1" = 'W_M'))+
  scale_y_continuous(limits = c(-0.21, 0.21),
                     breaks = seq(-0.2,0.2,0.1),
                     labels = c(-0.2,-0.1,0,0.1,0.2))+
  geom_hline(yintercept=0, lty=1)+
  ylab(NULL)+
  ylab(expression("Interactive effect on " * Delta * C[veg] * " (kg C " * m^-2 * " " * yr^-1 * ")"))+
  annotate('text', x=1, y=0.18, label='(e)', size=10)+
  theme_linedraw() +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.margin = unit(c(0.7,0.2,0,0.2), 'cm'),
        plot.title = element_text(size=24, face="bold",
                                  hjust = 0.5),
        # axis.text.x=element_text(size=20, angle = 35, vjust = 1.0, hjust = 1.0),
        axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=22),
        # axis.title.y=element_text(size=22, vjust = 1),
        axis.title.y = element_blank(),
        axis.ticks.length = unit(0.2, 'cm'),
        legend.text=element_text(size=23),
        legend.title=element_text(size=20),
        legend.position='none',
        legend.box = 'horizontal',
        legend.box.just = 'left',
        legend.key.size = unit(1.3, 'cm'))

 # NEP
 pi6 <-
  ggplot(data = subset(PE_effect_plot, Period == 'futu' & Effect == 'MDP_x_CO2' & Value == 'magnitude'))+
  geom_bar(stat = "identity",
           mapping=aes(ModName, NEP/1000, fill=ModName),
           position=position_dodge(preserve = 'single'),
           col="black",
           width = 0.6) +
  geom_errorbar(mapping=aes(x = ModName, 
                            ymax = NEP/1000 + NEP_sd/1000,
                            ymin = NEP/1000 - NEP_sd/1000,
                            col= 'black'),
                position=position_dodge(preserve = 'single'),
                col="black",
                width = 0.4)+
  geom_vline(xintercept=c(2.5), lty=3, linewidth=1)+
  scale_fill_manual(name="",
                    values=c('grey30', 'red', model.cols),
                    labels=c('M_M', 'M_M1', model.labels))+
  scale_x_discrete(name="",
                   # guide = guide_axis(n.dodge = 2),
                   labels=c("A_GDAYP" = "GDAYP",
                            "B_ELMV1" = "ELMV1",
                            "C_CABLP" = "CABLP",
                            "D_LPJGP" = "LPJGP",
                            "E_OCHDP" = "OCHDP",
                            "F_QUINC" = "QUINC",
                            "G_OCHDX" = "OCHDX",
                            "H_QUJSM" = "QUJSM",
                            "M_M" = "M_M",
                            "M_M1" = 'W_M'))+
  scale_y_continuous(limits = c(-0.21, 0.21),
                     breaks = seq(-0.2,0.2,0.1))+
  geom_hline(yintercept=0, lty=1)+
  ylab(NULL)+
  ylab(expression("Interactive effect on " * NEP * " (kg C " * m^-2 * " " * yr^-1 * ")"))+
  annotate('text', x=1, y=0.18, label='(f)', size=10)+
  theme_linedraw() +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.margin = unit(c(0.7,0.2,0.2,0.2), 'cm'),
        plot.title = element_text(size=24, face="bold",
                                  hjust = 0.5),
        axis.text.x=element_text(size=20, angle = 35, vjust = 1.0, hjust = 1.0),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=22),
        # axis.title.y=element_text(size=22, vjust = 1),
        axis.title.y = element_blank(),
        axis.ticks.length = unit(0.2, 'cm'),
        legend.text=element_text(size=23),
        legend.title=element_text(size=20),
        legend.position='none',
        legend.box = 'horizontal',
        legend.box.just = 'left',
        legend.key.size = unit(1.3, 'cm'))
  
  
  # MDP, 2023-2029
  # ----------------------------------------------------
  # GPP
  pi7 <-
  ggplot(data = subset(PE_effect_plot, Period == 'futu' & Effect == 'HIP_x_CO2' & Value == 'percentage'))+
    geom_bar(stat = "identity",
             mapping=aes(ModName, GPP, fill=ModName),
             position=position_dodge(preserve = 'single'),
             col="black",
             width = 0.6) +
    geom_errorbar(mapping=aes(x = ModName, 
                              ymax = GPP + GPP_sd,
                              ymin = GPP - GPP_sd,
                              col= 'black'),
                  position=position_dodge(preserve = 'single'),
                  col="black",
                  width = 0.4)+
    geom_vline(xintercept=c(2.5), lty=3, linewidth=1)+
    scale_fill_manual(name="",
                      values=c('grey30', 'red', model.cols),
                      labels=c('M_M', 'M_M1', model.labels))+
    scale_x_discrete(name="",
                     # guide = guide_axis(n.dodge = 2),
                     labels=c("A_GDAYP" = "GDAYP",
                              "B_ELMV1" = "ELMV1",
                              "C_CABLP" = "CABLP",
                              "D_LPJGP" = "LPJGP",
                              "E_OCHDP" = "OCHDP",
                              "F_QUINC" = "QUINC",
                              "G_OCHDX" = "OCHDX",
                              "H_QUJSM" = "QUJSM",
                              "M_M" = "M_M",
                              "M_M1" = 'W_M'))+
    scale_y_continuous(limits = c(-20, 10),
                       breaks = seq(-20,10,10))+
    geom_hline(yintercept=0, lty=1)+
    ylab(NULL)+
    ylab(expression("Interactive effect on " * GPP * " (%)"))+
    ggtitle('HP, 2023-2029') +
    annotate('text', x=1, y=9, label='(g)', size=10)+
    theme_linedraw() +
    theme(panel.grid.minor=element_blank(),
            panel.grid.major=element_blank(),
            plot.margin = unit(c(1,1.5,0,0.2), 'cm'),
            plot.title = element_text(size=30, 
                                      hjust = 0.5),
            # axis.text.x=element_text(size=20, angle = 35, vjust = 1.0, hjust = 1.0),
            axis.text.x = element_blank(),
            axis.title.x=element_blank(),
            axis.text.y=element_text(size=22),
            # axis.title.y=element_text(size=22, vjust = 1),
            axis.title.y = element_blank(),
            axis.ticks.length = unit(0.2, 'cm'),
            legend.text=element_text(size=23),
            legend.title=element_text(size=20),
            legend.position='none',
            legend.box = 'horizontal',
            legend.box.just = 'left',
            legend.key.size = unit(1.3, 'cm'))
  
  # biomass growth
  pi8 <-
  ggplot(data = subset(PE_effect_plot, Period == 'futu' & Effect == 'HIP_x_CO2' & Value == 'magnitude'))+
    geom_bar(stat = "identity",
             mapping=aes(ModName, deltaCveg/1000, fill=ModName),
             position=position_dodge(preserve = 'single'),
             col="black",
             width = 0.6) +
    geom_errorbar(mapping=aes(x = ModName, 
                              ymax = deltaCveg/1000 + deltaCveg_sd/1000,
                              ymin = deltaCveg/1000 - deltaCveg_sd/1000,
                              col= 'black'),
                  position=position_dodge(preserve = 'single'),
                  col="black",
                  width = 0.4)+
    geom_vline(xintercept=c(2.5), lty=3, linewidth=1)+
    scale_fill_manual(name="",
                      values=c('grey30', 'red', model.cols),
                      labels=c('M_M', 'M_M1', model.labels))+
    scale_x_discrete(name="",
                     # guide = guide_axis(n.dodge = 2),
                     labels=c("A_GDAYP" = "GDAYP",
                              "B_ELMV1" = "ELMV1",
                              "C_CABLP" = "CABLP",
                              "D_LPJGP" = "LPJGP",
                              "E_OCHDP" = "OCHDP",
                              "F_QUINC" = "QUINC",
                              "G_OCHDX" = "OCHDX",
                              "H_QUJSM" = "QUJSM",
                              "M_M" = "M_M",
                              "M_M1" = 'W_M'))+
    scale_y_continuous(limits = c(-0.21, 0.21),
                       breaks = seq(-0.2,0.2,0.1),
                       labels = c(-0.2,-0.1,0,0.1,0.2))+
    geom_hline(yintercept=0, lty=1)+
    ylab(NULL)+
    ylab(expression("Interactive effect on " * Delta * C[veg] * " (kg C " * m^-2 * " " * yr^-1 * ")"))+
    annotate('text', x=1, y=0.18, label='(h)', size=10)+
    theme_linedraw() +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          plot.margin = unit(c(0.7,1.5,0,0.2), 'cm'),
          plot.title = element_text(size=24, face="bold",
                                    hjust = 0.5),
          # axis.text.x=element_text(size=20, angle = 35, vjust = 1.0, hjust = 1.0),
          axis.text.x = element_blank(),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=22),
          # axis.title.y=element_text(size=22, vjust = 1),
          axis.title.y = element_blank(),
          axis.ticks.length = unit(0.2, 'cm'),
          legend.text=element_text(size=23),
          legend.title=element_text(size=20),
          legend.position='none',
          legend.box = 'horizontal',
          legend.box.just = 'left',
          legend.key.size = unit(1.3, 'cm'))
  
  # NEP
  pi9 <-
  ggplot(data = subset(PE_effect_plot, Period == 'futu' & Effect == 'HIP_x_CO2' & Value == 'magnitude'))+
    geom_bar(stat = "identity",
             mapping=aes(ModName, NEP/1000, fill=ModName),
             position=position_dodge(preserve = 'single'),
             col="black",
             width = 0.6) +
    geom_errorbar(mapping=aes(x = ModName, 
                              ymax = NEP/1000 + NEP_sd/1000,
                              ymin = NEP/1000 - NEP_sd/1000,
                              col= 'black'),
                  position=position_dodge(preserve = 'single'),
                  col="black",
                  width = 0.4)+
    geom_vline(xintercept=c(2.5), lty=3, linewidth=1)+
    scale_fill_manual(name="",
                      values=c('grey30', 'red', model.cols),
                      labels=c('M_M', 'M_M1', model.labels))+
    scale_x_discrete(name="",
                     # guide = guide_axis(n.dodge = 2),
                     labels=c("A_GDAYP" = "GDAYP",
                              "B_ELMV1" = "ELMV1",
                              "C_CABLP" = "CABLP",
                              "D_LPJGP" = "LPJGP",
                              "E_OCHDP" = "OCHDP",
                              "F_QUINC" = "QUINC",
                              "G_OCHDX" = "OCHDX",
                              "H_QUJSM" = "QUJSM",
                              "M_M" = "M_M",
                              "M_M1" = 'W_M'))+
    scale_y_continuous(limits = c(-0.21, 0.21),
                       breaks = seq(-0.2,0.2,0.12))+
    geom_hline(yintercept=0, lty=1)+
    ylab(NULL)+
    ylab(expression("Interactive effect on " * NEP * " (kg C " * m^-2 * " " * yr^-1 * ")"))+
    annotate('text', x=1, y=0.18, label='(i)', size=10)+
    theme_linedraw() +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          plot.margin = unit(c(0.7,1.5,0.2,0.2), 'cm'),
          plot.title = element_text(size=24, face="bold",
                                    hjust = 0.5),
          axis.text.x=element_text(size=20, angle = 35, vjust = 1.0, hjust = 1.0),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=22),
          # axis.title.y=element_text(size=22, vjust = 1),
          axis.title.y = element_blank(),
          axis.ticks.length = unit(0.2, 'cm'),
          legend.text=element_text(size=23),
          legend.title=element_text(size=20),
          legend.position='none',
          legend.box = 'horizontal',
          legend.box.just = 'left',
          legend.key.size = unit(1.3, 'cm'))
  
  
  plot_interaction <- 
    pi1 + pi2 + pi3 + 
    pi4 + pi5 + pi6 + 
    pi7 + pi8 + pi9 + 
    plot_layout(ncol = 3, byrow = F)
}




#############################################################################
#### Modelling CO2 x P at EucFACE: multi-model intercomparison           ####
#### P fertilization effects                                             ####
#### Created by Bin Wang, 2025-04-27                                     ####
#### Email: wangbin1992@zju.edu.cn                                       ####
#############################################################################




### clear space
rm(list = ls());



#### preparations
#============================================
##packages
library(doBy)
library(ggplot2)
library(ggpattern)
library(cowplot)
library(gridExtra)
library(patchwork)
library(reshape2)
library(ggplotify)
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




#### data prepare
# ===================================================================================

## read annual datasets
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

## calculate CUE, Cveg, deltaCveg, canopy P, leaf P concentration
ambDF$CUE <- ambDF$NPP/ambDF$GPP
eleDF$CUE <- eleDF$NPP/eleDF$GPP
ambDF$Cveg <- rowSums(data.frame(ambDF$CL,
                                 ambDF$CW,
                                 ambDF$CFR,
                                 ambDF$CCR,
                                 ambDF$CSTOR),
                      na.rm = T)
eleDF$Cveg <- rowSums(data.frame(eleDF$CL,
                                 eleDF$CW,
                                 eleDF$CFR,
                                 eleDF$CCR,
                                 eleDF$CSTOR),
                      na.rm = T)
ambDF$deltaCveg <- rowSums(data.frame(ambDF$deltaCL,
                                      ambDF$deltaCW,
                                      ambDF$deltaCFR,
                                      ambDF$deltaCCR,
                                      ambDF$deltaCSTOR),
                           na.rm = TRUE)
eleDF$deltaCveg <- rowSums(data.frame(eleDF$deltaCL,
                                      eleDF$deltaCW,
                                      eleDF$deltaCFR,
                                      eleDF$deltaCCR,
                                      eleDF$deltaCSTOR),
                           na.rm = TRUE)
ambDF$PCANOPY <- ambDF$PL*ambDF$LAI # P mass of leaves multiply by LAI
eleDF$PCANOPY <- eleDF$PL*eleDF$LAI # P mass of leaves multiply by LAI
ambDF$PLCON <- ambDF$PL/ambDF$LAI # P mass of leaves divided by LAI
eleDF$PLCON <- eleDF$PL/eleDF$LAI # P mass of leaves divided by LAI

## organize to make neat
ambDF$CO2TRT <- 'ambient'
eleDF$CO2TRT <- 'elevated'
ambDF <- ambDF[c(1,2,157,151,3:150,152:156)]
eleDF <- eleDF[c(1,2,157,151,3:150,152:156)]





#### data processing and plot, HP, 2020-2022
# ===================================================================

## data processing
# =================================
{
  ## papare HP, CO2, and HP + CO2 effects
  NOP_amb <- subset(ambDF, (YEAR >= 2020 & YEAR <= 2022) &
                      PTRT == 'NOP',
                    select = c('ModName', 'YEAR', 'GPP', 'deltaCveg', 'NEP', 
                               'CUE', 'NPP', 'deltaCL', 'deltaCW', 'deltaCFR', 'deltaCCR', 'deltaCSTOR', 'deltaCSOIL'))
  HIP_amb <- subset(ambDF, (YEAR >= 2020 & YEAR <= 2022) &
                      PTRT == 'HIP',
                    select = c('ModName', 'YEAR', 'GPP', 'deltaCveg', 'NEP', 
                               'CUE', 'NPP', 'deltaCL', 'deltaCW', 'deltaCFR', 'deltaCCR', 'deltaCSTOR', 'deltaCSOIL'))
  NOP_ele <- subset(eleDF, (YEAR >= 2020 & YEAR <= 2022) &
                      PTRT == 'NOP',
                    select = c('ModName', 'YEAR', 'GPP', 'deltaCveg', 'NEP', 
                               'CUE', 'NPP', 'deltaCL', 'deltaCW', 'deltaCFR', 'deltaCCR', 'deltaCSTOR', 'deltaCSOIL'))
  HIP_ele <- subset(eleDF, (YEAR >= 2020 & YEAR <= 2022) &
                      PTRT == 'HIP',
                    select = c('ModName', 'YEAR', 'GPP', 'deltaCveg', 'NEP', 
                               'CUE', 'NPP', 'deltaCL', 'deltaCW', 'deltaCFR', 'deltaCCR', 'deltaCSTOR', 'deltaCSOIL'))
  
  # calculate P effect, magnitude
  d1 <- dim(NOP_amb)[2]
  Peff_amb <- NOP_amb # P effect
  Peff_amb[, 3:d1] <- HIP_amb[, 3:d1] - NOP_amb[, 3:d1]
  Peff_ele <- NOP_amb # P effect
  Peff_ele[, 3:d1] <- HIP_ele[, 3:d1] - NOP_ele[, 3:d1]
  
  # calculate P effect, %
  Peff_amb1 <- NOP_amb # P effect
  Peff_amb1[, 3:d1] <- (HIP_amb[, 3:d1] - NOP_amb[, 3:d1])/NOP_amb[, 3:d1]*100
  Peff_ele1 <- NOP_amb # P effect
  Peff_ele1[, 3:d1] <- (HIP_ele[, 3:d1] - NOP_ele[, 3:d1])/NOP_ele[, 3:d1]*100
  
  Peff_amb$CO2TRT <- Peff_amb1$CO2TRT <- 'amb'
  Peff_ele$CO2TRT <- Peff_ele1$CO2TRT <- 'ele'
  Peff_amb$SIZE <- Peff_ele$SIZE <- 'mag'
  Peff_amb1$SIZE <- Peff_ele1$SIZE <- 'per'
  
  PE_effect <- rbind(Peff_amb, rbind(Peff_amb1, rbind(Peff_ele, Peff_ele1)))
  PE_effect_model <- summaryBy(.~ModName+CO2TRT+SIZE, data = PE_effect[-2], FUN=mean, keep.names = T, na.rm=TRUE)
}



## plot
# =================================
{
  
  # prepare data
  PE_effect_mip1 <- summaryBy(.~CO2TRT+SIZE, data = PE_effect_model[2:14], FUN=mean, keep.names = T, na.rm=T)
  PE_effect_mip11 <- summaryBy(.~CO2TRT+SIZE, data = PE_effect_model[2:14], FUN=sd, keep.names = T, na.rm=T)
  
  PE_effect_mip_1 <- melt(PE_effect_mip1, id=c('CO2TRT', 'SIZE'))
  names(PE_effect_mip_1)[4] <- 'mean'
  PE_effect_mip_11 <- melt(PE_effect_mip11, id=c('CO2TRT', 'SIZE'))
  names(PE_effect_mip_11)[4] <- 'sd'
  
  PE_effect_mip <- merge(PE_effect_mip_1, PE_effect_mip_11)
  PE_effect_mip$variable <- factor(PE_effect_mip$variable,
                                   levels = c('GPP', 'NPP', 'CUE', 
                                              'deltaCveg', 'deltaCL', 'deltaCW', 'deltaCFR', 'deltaCCR', 'deltaCSTOR','deltaCSOIL', 'NEP'))
  
  PE_effect_model1 <- melt(PE_effect_model, id=c('ModName', 'CO2TRT', 'SIZE'))
  PE_effect_model1$variable <- factor(PE_effect_model1$variable,
                                      levels = c('GPP', 'NPP', 'CUE', 
                                                 'deltaCveg', 'deltaCL', 'deltaCW', 'deltaCFR', 'deltaCCR', 'deltaCSTOR', 'deltaCSOIL', 'NEP'))
  
  # percentage
  p1 <-
    ggplot()+
    geom_bar(stat = "identity",
             data = subset(PE_effect_mip, SIZE == 'per' & 
                             (variable == 'GPP' |
                                variable == 'NPP' | 
                                variable == 'CUE')),
             mapping=aes(x=variable, y=mean, group=CO2TRT, fill=CO2TRT),
             position=position_dodge(preserve = 'single'),
             col="black",
             width = 0.75,
             alpha=0.5)+
    geom_errorbar(data = subset(PE_effect_mip, SIZE == 'per' & 
                                  (variable == 'GPP' |
                                     variable == 'NPP' | 
                                     variable == 'CUE')),
                  mapping=aes(x = variable, 
                              ymax = mean + sd,
                              ymin = mean - sd,
                              group=CO2TRT,
                              col= 'black'),
                  position=position_dodge(width = 0.75),
                  col="black",
                  width = 0.38,
                  linewidth=1.5)+
    geom_point(data = subset(PE_effect_model1, SIZE == 'per' & 
                               (variable == 'GPP' |
                                  variable == 'NPP' | 
                                  variable == 'CUE')),
               aes(x=variable, y=value, group=CO2TRT, fill=ModName),
               position=position_dodge(width = 0.75),
               pch=21,
               size=6,
               alpha=0.75)+
    scale_x_discrete(name='',
                     labels = c('GPP' = 'GPP',
                                'NPP' = 'NPP',
                                'CUE' = 'CUE'))+
    scale_y_continuous(limits = c(-45,125),
                       breaks = seq(-40,120,40))+
    scale_fill_manual(values = c(model.cols[1],
                                 'blue',
                                 model.cols[2:5],
                                 'red',
                                 model.cols[6:8]),
                      labels = c(model.labels[1],
                                 'amb' = expression(aCO[2]),
                                 model.labels[2:5],
                                 'ele' = expression(eCO[2]),
                                 model.labels[6:8]))+
    geom_hline(yintercept = 0, lty=1)+
    ylab(expression("P effect (%)"  ))+
    annotate('text', x=0.70, y=120, label='(a)', size=11)+
    theme_linedraw() +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          plot.margin = unit(c(1,1,4,0.4), 'cm'),
          axis.ticks.length = unit(0.2, 'cm'),
          plot.title = element_blank(),
          axis.text.x=element_text(size=24, vjust = 0),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=24),
          axis.title.y=element_text(size=24, vjust = 1),
          legend.text=element_text(size=24),
          # legend.position=c(0.8,0.85),
          legend.position = 'none',
          legend.box = 'horizontal',
          legend.box.just = 'left',
          legend.spacing.x = unit(0.6, 'cm'),
          legend.key.width = unit(1.4, 'cm'),
          legend.key.height = unit(1.2, 'cm'))
  
  # magnitude
  p2 <-
    ggplot()+
    geom_bar(stat = "identity",
             data = subset(PE_effect_mip, SIZE == 'mag' & 
                             (variable == 'deltaCveg' |
                                variable == 'deltaCL' |
                                variable == 'deltaCW' |
                                variable == 'deltaCFR' |
                                variable == 'deltaCCR' |
                                variable == 'deltaCSTOR' |
                                variable == 'deltaCSOIL' |
                                variable == 'NEP' )),
             mapping=aes(x=variable, y=mean/1000, group=CO2TRT, fill=CO2TRT),
             position=position_dodge(preserve = 'single'),
             col="black",
             width = 0.75,
             alpha=0.5)+
    geom_errorbar(data = subset(PE_effect_mip, SIZE == 'mag' & 
                                  (variable == 'deltaCveg' |
                                     variable == 'deltaCL' |
                                     variable == 'deltaCW' |
                                     variable == 'deltaCFR' |
                                     variable == 'deltaCCR' |
                                     variable == 'deltaCSTOR' |
                                     variable == 'deltaCSOIL' |
                                     variable == 'NEP' )),
                  mapping=aes(x = variable, 
                              ymax = mean/1000 + sd/1000,
                              ymin = mean/1000 - sd/1000,
                              group=CO2TRT,
                              col= 'black'),
                  position=position_dodge(width = 0.75),
                  col="black",
                  width = 0.38,
                  linewidth=1.5)+
    geom_point(data = subset(PE_effect_model1, SIZE == 'mag' & 
                               (variable == 'deltaCveg' |
                                  variable == 'deltaCL' |
                                  variable == 'deltaCW' |
                                  variable == 'deltaCFR' |
                                  variable == 'deltaCCR' |
                                  variable == 'deltaCSTOR' |
                                  variable == 'deltaCSOIL' |
                                  variable == 'NEP' )),
               aes(x=variable, y=value/1000, group=CO2TRT, fill=ModName),
               position=position_dodge(width = 0.75),
               pch=21,
               size=6,
               alpha=0.75)+
    # scale_y_continuous(limits = c(-0.22,0.6),
    #                    breaks = seq(-0.2,0.6,0.2))+
    scale_x_discrete(name='',
                     labels = c('deltaCveg' = expression(Delta * C[veg]),
                                'deltaCL' = expression(Delta * C[leaf]),
                                'deltaCW' = expression(Delta * C[wood]),
                                'deltaCFR' = expression(Delta * C[froot]),
                                'deltaCCR' = expression(Delta * C[croot]),
                                'deltaCSTOR' = expression(Delta * C[store]),
                                'deltaCSOIL' = expression(Delta * C[soil]),
                                'NEP' = 'NEP'))+
    scale_fill_manual(values = c(model.cols[1],
                                 'blue',
                                 model.cols[2:5],
                                 'red',
                                 model.cols[6:8]),
                      labels = c(model.labels[1],
                                 'amb' = expression(aCO[2]),
                                 model.labels[2:5],
                                 'ele' = expression(eCO[2]),
                                 model.labels[6:8]))+
    geom_hline(yintercept = 0, lty=1)+
    # ylab('P effect')+
    ylab(expression("P effect (kg C " * m^-2 * " " * yr^-1 * ")"  ))+
    annotate('text', x=0.70, y=0.51, label='(b)', size=11)+
    theme_linedraw() +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          plot.margin = unit(c(1,1,4,1), 'cm'),
          axis.ticks.length = unit(0.2, 'cm'),
          plot.title = element_blank(),
          axis.text.x=element_text(size=24, vjust = 0),
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=24),
          axis.title.y=element_text(size=24, vjust = 1),
          legend.text=element_text(size=24),
          legend.title=element_blank(),
          legend.position='none',
          legend.box = 'horizontal',
          legend.box.just = 'left',
          legend.key.size = unit(1.8, 'cm'))

}


## combine
(p1 + p2) + 
  plot_layout(ncol = 2,
              widths = c(1.5,4.5))


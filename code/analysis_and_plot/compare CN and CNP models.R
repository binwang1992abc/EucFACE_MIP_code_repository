

#############################################################################
#### Modelling CO2 x P at EucFACE: multi-model intercomparison           ####
#### Comparison of CN and CNP models, GDAY and LPJ-GUESS                 ####
#### Created by Bin Wang, 2025-04-27                                     ####
#### Email: wangbin1992@zju.edu.cn                                       ####
#############################################################################



rm(list = ls());




#### preparation
# ======================================================

## packages
library(doBy)
library(ggplot2)
library(ggpattern)
library(cowplot)
library(gridExtra)
library(patchwork)
library(reshape2)
library(RColorBrewer)

## define climate
climate.scenario <- 'FIX' # fixed climate in wet

## set pathways of input and output
setwd('D:/Analysis/EucFACE_MIP/EucFACE_MIP_P_x_CO2/')

## color set
GreensPalette <- rev(brewer.pal(n = 9, name = "Greens"))
YlOrRdPalette <- rev(brewer.pal(n = 9, name = "YlOrRd"))
SpectralPalette <- brewer.pal(n = 5, name = "Spectral")
Diverge_hsv_Palette <- colorspace::diverge_hcl(5)




#### data preparation, ambient dataset
# ========================================================================

## read annual datasets, ambient
#-------------------------------------------------------------------------
ambNOP <- readRDS(paste0(getwd(),
                         "/data/compile_output/MIP_ALL_", 
                         climate.scenario, "_NOP_AMB_annual.rds"))
ambMDP <- readRDS(paste0(getwd(),
                         "/data/compile_output/MIP_ALL_", 
                         climate.scenario, "_MDP_AMB_annual.rds"))
ambHIP <- readRDS(paste0(getwd(),
                         "/data/compile_output/MIP_ALL_", 
                         climate.scenario, "_HIP_AMB_annual.rds"))

## add P treatment, combine
ambNOP$PTRT <- "NOP"
ambMDP$PTRT <- "MDP"
ambHIP$PTRT <- "HIP"
ambDF <- rbind(ambNOP, rbind(ambMDP, ambHIP))

## ignore climate variables (fixed variable)
ambDF[,c("PAR", "TAIR", "TSOIL", "VPD", "PREC")] <- NULL

## ignore NAs
ambDF[ambDF<=-999] <- NA

## select N-contained model
ambDF <- subset(ambDF, (ModName == 'I_GDAYN' | ModName == 'J_LPJGN' | ModName == 'A_GDAYP' | ModName == 'D_LPJGP'))

## remove 2012
ambDF <- subset(ambDF, YEAR >= 2013)

## calculate biomass C and its growth
ambDF$Cveg <- rowSums(data.frame(ambDF$CL, ambDF$CW, ambDF$CFR, ambDF$CCR, ambDF$CSTOR), na.rm = T)
ambDF$deltaCveg <- rowSums(data.frame(ambDF$deltaCL, ambDF$deltaCW, ambDF$deltaCFR,
                                      ambDF$deltaCCR, ambDF$deltaCSTOR), na.rm = TRUE)


ambDF_PA <- subset(ambDF, YEAR >= 2020 & YEAR <=2022)
ambDF_PA <- ambDF_PA[-which(ambDF_PA$ModName == 'I_GDAYN' & ambDF_PA$PTRT == 'MDP'), ]
ambDF_PA <- ambDF_PA[-which(ambDF_PA$ModName == 'I_GDAYN' & ambDF_PA$PTRT == 'HIP'), ]
ambDF_PA <- ambDF_PA[-which(ambDF_PA$ModName == 'J_LPJGN' & ambDF_PA$PTRT == 'MDP'), ]
ambDF_PA <- ambDF_PA[-which(ambDF_PA$ModName == 'J_LPJGN' & ambDF_PA$PTRT == 'HIP'), ]

ambDF_PA1 <- ambDF_PA
ambDF_PA1$ModPTRT <- factor(paste(ambDF_PA1$ModName, ambDF_PA1$PTRT, sep = '_'))
CN_ANA <- summaryBy(.~PTRT+ModName+ModPTRT, data = ambDF_PA1[-2], FUN = c(mean, sd), keep.names = T)



## read annual datasets, elevated
#-----------------------------------------------------------------------------------
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
eleNOP$PTRT <- "NOP"
eleMDP$PTRT <- "MDP"
eleHIP$PTRT <- "HIP"
eleDF <- rbind(eleNOP, rbind(eleMDP, eleHIP))

## ignore climate variables (fixed variable)
eleDF[,c("PAR", "TAIR", "TSOIL", "VPD", "PREC")] <- NULL

## ignore NAs
eleDF[eleDF<=-999] <- NA

## select N-contained model
eleDF <- subset(eleDF, (ModName == 'I_GDAYN' | ModName == 'J_LPJGN' | ModName == 'A_GDAYP' | ModName == 'D_LPJGP'))

## remove 2012
eleDF <- subset(eleDF, YEAR >= 2013)

## calculate biomass C and its growth
eleDF$Cveg <- rowSums(data.frame(eleDF$CL, eleDF$CW, eleDF$CFR, eleDF$CCR, eleDF$CSTOR), na.rm = T)
eleDF$deltaCveg <- rowSums(data.frame(eleDF$deltaCL, eleDF$deltaCW, eleDF$deltaCFR,
                                      eleDF$deltaCCR, eleDF$deltaCSTOR), na.rm = TRUE)


eleDF_PA <- subset(eleDF, YEAR >= 2020 & YEAR <=2022)
eleDF_PA <- eleDF_PA[-which(eleDF_PA$ModName == 'I_GDAYN' & eleDF_PA$PTRT == 'MDP'), ]
eleDF_PA <- eleDF_PA[-which(eleDF_PA$ModName == 'I_GDAYN' & eleDF_PA$PTRT == 'HIP'), ]
eleDF_PA <- eleDF_PA[-which(eleDF_PA$ModName == 'J_LPJGN' & eleDF_PA$PTRT == 'MDP'), ]
eleDF_PA <- eleDF_PA[-which(eleDF_PA$ModName == 'J_LPJGN' & eleDF_PA$PTRT == 'HIP'), ]



## calculate CO2 effect
#-----------------------------------------------------------------------------------

## NOP
amb_NOP <- subset(ambDF_PA, PTRT == 'NOP', select = c('ModName', 'YEAR', 'PTRT', 'GPP', 'deltaCveg', 'NEP'))
ele_NOP <- subset(eleDF_PA, PTRT == 'NOP', select = c('ModName', 'YEAR', 'PTRT', 'GPP', 'deltaCveg', 'NEP'))
CO2Dif_NOP <- amb_NOP
CO2Dif_NOP[4:6] <- ele_NOP[4:6] - amb_NOP[4:6]

## MDP and HIP
ele_MDP <- subset(eleDF_PA, PTRT == 'MDP', select = c('ModName', 'YEAR', 'PTRT', 'GPP', 'deltaCveg', 'NEP'))
ele_HIP <- subset(eleDF_PA, PTRT == 'HIP', select = c('ModName', 'YEAR', 'PTRT', 'GPP', 'deltaCveg', 'NEP'))
amb_NOP_1 <- subset(ambDF_PA, PTRT == 'NOP' & !(ModName == 'I_GDAYN' | ModName == 'J_LPJGN'), select = c('ModName', 'YEAR', 'PTRT', 'GPP', 'deltaCveg', 'NEP'))

CO2Dif_MDP <- ele_MDP
CO2Dif_MDP[4:6] <- ele_MDP[4:6] - amb_NOP_1[4:6]

CO2Dif_HIP <- ele_HIP
CO2Dif_HIP[4:6] <- ele_HIP[4:6] - amb_NOP_1[4:6]

CO2Dif <- rbind(CO2Dif_NOP, rbind(CO2Dif_MDP, CO2Dif_HIP))
CO2Dif$ModPTRT <- factor(paste(CO2Dif$ModName, CO2Dif$PTRT, sep = '_'))

CN_ANA1 <- summaryBy(.~PTRT+ModName+ModPTRT, data = CO2Dif[-2], FUN = c(mean, sd), keep.names = T)



### plot
# ============================================================================

## comparation of simualted GPP, biomass growth, and NEP
# ---------------------------------------------------------------------------------------
CN_ANA_GDAY <- subset(CN_ANA, ModName == 'A_GDAYP' | ModName == 'I_GDAYN')
CN_ANA_GDAY$ModPTRT <- factor(CN_ANA_GDAY$ModPTRT, 
                              levels = c('I_GDAYN_NOP', 'A_GDAYP_NOP', 'A_GDAYP_MDP', 'A_GDAYP_HIP'),
                              labels = c('GDAYN', 'GDAYP_NOP', 'GDAYP_MDP', 'GDAYP_HIP'))
p1 <-
  ggplot()+
  geom_bar(stat = "identity", 
           data = CN_ANA_GDAY,
           mapping=aes(ModPTRT, GPP.mean/1000, fill=ModPTRT),
           position=position_dodge(preserve = 'single'), 
           alpha = 0.7,
           col="black",
           width = 0.5)+
  geom_errorbar(stat = "identity",
                data = CN_ANA_GDAY,
                mapping=aes(ModPTRT,
                            ymax=GPP.mean/1000+GPP.sd/1000,
                            ymin=GPP.mean/1000-GPP.sd/1000,
                            group=PTRT),
                # position=position_dodge(preserve = 'single', width=0.7),
                linewidth=1,
                width=0.3)+
  scale_fill_manual(name="",
                    values=c("GDAYN"="cyan4",
                             "GDAYP_NOP"="grey70",
                             "GDAYP_MDP"="goldenrod1",
                             "GDAYP_HIP"="seagreen"))+
  scale_x_discrete(name='',
                   labels = c("GDAYN"="CN",
                              "GDAYP_NOP"="aP",
                              "GDAYP_MDP"="MP",
                              "GDAYP_HIP"="HP"))+
  scale_y_continuous(limits = c(0, 3),
                     breaks = seq(0, 3, 0.5))+
  geom_hline(yintercept=0, lty=1)+
  ylab(expression(GPP * "  ( kg C " * m^-2 * " " * yr^-1 * " )"))+
  ggtitle('GDAY')+
  annotate('text', x=0.73, y=2.85, label='(a)', size=10)+
  theme_linedraw() +
  guides(fill = 'none')+
  theme(text = element_text(family = 'serif'),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.margin = unit(c(0.5,0.8,0,0.5), 'cm'),
        plot.title = element_text(size=32,
                                  hjust = 0.5),
        # axis.text.x=element_text(size=22, vjust = -0.5),
        axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=22),
        axis.title.y=element_text(size=24, vjust = 2.5,
                                  margin = unit(c(0.28,0.28,0.28,0.28), 'cm')),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.position=c(0.2,0.85),
        legend.box = 'horizontal',
        legend.box.just = 'left',
        legend.key.size = unit(0.8, 'cm'))

p2 <-
  ggplot()+
  geom_bar(stat = "identity", 
           data = CN_ANA_GDAY,
           mapping=aes(ModPTRT, deltaCveg.mean/1000, fill=ModPTRT),
           position=position_dodge(preserve = 'single'), 
           alpha = 0.7,
           col="black",
           width = 0.5)+
  geom_errorbar(stat = "identity",
                data = CN_ANA_GDAY,
                mapping=aes(ModPTRT,
                            ymax=deltaCveg.mean/1000+deltaCveg.sd/1000,
                            ymin=deltaCveg.mean/1000-deltaCveg.sd/1000,
                            group=PTRT),
                # position=position_dodge(preserve = 'single', width=0.7),
                linewidth=1,
                width=0.3)+
  scale_fill_manual(name="",
                    values=c("GDAYN"="cyan4",
                             "GDAYP_NOP"="grey70",
                             "GDAYP_MDP"="goldenrod1",
                             "GDAYP_HIP"="seagreen"))+
  scale_x_discrete(name='',
                   labels = c("GDAYN"="CN",
                              "GDAYP_NOP"="aP",
                              "GDAYP_MDP"="MP",
                              "GDAYP_HIP"="HP"))+
  scale_y_continuous(limits = c(-0.05, 0.4),
                     breaks = seq(0, 0.4, 0.1))+
  geom_hline(yintercept=0, lty=1)+
  ylab(expression(Delta * C[veg] * "  ( kg C " * m^-2 * " " * yr^-1 * " )"))+
  ggtitle('')+
  annotate('text', x=0.73, y=0.38, label='(e)', size=10)+
  theme_linedraw() +
  guides(fill = FALSE)+
  theme(text = element_text(family = 'serif'),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.margin = unit(c(0,0.8,0,0.5), 'cm'),
        plot.title = element_text(size=16, face="bold",
                                  hjust = 0.5),
        # axis.text.x=element_text(size=22, vjust = -0.5),
        axis.text.x=element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=22),
        axis.title.y=element_text(size=24, vjust = 2.5,
                                  margin = unit(c(0.28,0.28,0.28,0.28), 'cm')),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.position='none',
        legend.box = 'horizontal',
        legend.box.just = 'left',
        legend.key.size = unit(0.8, 'cm'))

p3 <-
  ggplot()+
  geom_bar(stat = "identity", 
           data = CN_ANA_GDAY,
           mapping=aes(ModPTRT, NEP.mean/1000, fill=ModPTRT),
           position=position_dodge(preserve = 'single'), 
           alpha = 0.7,
           col="black",
           width = 0.5)+
  geom_errorbar(stat = "identity",
                data = CN_ANA_GDAY,
                mapping=aes(ModPTRT,
                            ymax=NEP.mean/1000+NEP.sd/1000,
                            ymin=NEP.mean/1000-NEP.sd/1000,
                            group=PTRT),
                # position=position_dodge(preserve = 'single', width=0.7),
                linewidth=1,
                width=0.35)+
  scale_fill_manual(name="",
                    values=c("GDAYN"="cyan4",
                             "GDAYP_NOP"="grey70",
                             "GDAYP_MDP"="goldenrod1",
                             "GDAYP_HIP"="seagreen"))+
  scale_x_discrete(name='',
                   labels = c("GDAYN"="CN",
                              "GDAYP_NOP"="aP",
                              "GDAYP_MDP"="MP",
                              "GDAYP_HIP"="HP"))+
  scale_y_continuous(limits = c(-0.05, 0.5),
                     breaks = seq(0, 0.5, 0.1))+
  geom_hline(yintercept=0, lty=1)+
  ylab(expression(NEP * "  ( kg C " * m^-2 * " " * yr^-1 * " )"))+
  ggtitle('')+
  annotate('text', x=0.73, y=0.48, label='(i)', size=10)+
  theme_linedraw() +
  guides(fill = FALSE)+
  theme(text = element_text(family = 'serif'),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.margin = unit(c(0,0.8,0.5,0.5), 'cm'),
        plot.title = element_text(size=16, face="bold",
                                  hjust = 0.5),
        axis.text.x=element_text(size=28, vjust = -0.5),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=22),
        axis.title.y=element_text(size=24, vjust = 2.5,
                                  margin = unit(c(0.28,0.28,0.28,0.28), 'cm')),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.position='none',
        legend.box = 'horizontal',
        legend.box.just = 'left',
        legend.key.size = unit(0.8, 'cm'))


## LPJ-GUESS
CN_ANA_LPJ <- subset(CN_ANA, ModName == 'D_LPJGP' | ModName == 'J_LPJGN')
CN_ANA_LPJ$ModPTRT <- factor(CN_ANA_LPJ$ModPTRT, 
                             levels = c('J_LPJGN_NOP', 'D_LPJGP_NOP', 'D_LPJGP_MDP', 'D_LPJGP_HIP'),
                             labels = c('LPJGN', 'LPJGP_NOP', 'LPJGP_MDP', 'LPJGP_HIP'))
p4 <-
  ggplot()+
  geom_bar(stat = "identity", 
           data = CN_ANA_LPJ,
           mapping=aes(ModPTRT, GPP.mean/1000, fill=ModPTRT),
           position=position_dodge(preserve = 'single'), 
           alpha = 0.7,
           col="black",
           width = 0.5)+
  geom_errorbar(stat = "identity",
                data = CN_ANA_LPJ,
                mapping=aes(ModPTRT,
                            ymax=GPP.mean/1000+GPP.sd/1000,
                            ymin=GPP.mean/1000-GPP.sd/1000,
                            group=PTRT),
                # position=position_dodge(preserve = 'single', width=0.7),
                linewidth=1,
                width=0.3)+
  scale_fill_manual(name="",
                    values=c("LPJGN"="cyan4",
                             "LPJGP_NOP"="grey70",
                             "LPJGP_MDP"="goldenrod1",
                             "LPJGP_HIP"="seagreen"))+
  scale_x_discrete(name='',
                   labels = c("LPJGN"="CN",
                              "LPJGP_NOP"="aP",
                              "LPJGP_MDP"="MP",
                              "LPJGP_HIP"="HP"))+
  scale_y_continuous(limits = c(0, 2.5),
                     breaks = seq(0, 2.5, 0.5))+
  geom_hline(yintercept=0, lty=1)+
  ylab(expression(GPP * "  ( kg C " * m^-2 * " " * yr^-1 * " )"))+
  ggtitle('LPJ-GUESS')+
  annotate('text', x=0.73, y=2.4, label='(c)', size=10)+
  theme_linedraw() +
  guides(fill = FALSE)+
  theme(text = element_text(family = 'serif'),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.margin = unit(c(0.5,0.8,0,0.5), 'cm'),
        plot.title = element_text(size=34,
                                  hjust = 0.5),
        # axis.text.x=element_text(size=22, vjust = -0.5),
        axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=22),
        axis.title.y=element_text(size=24, vjust = 2.5,
                                  margin = unit(c(0.28,0.28,0.28,0.28), 'cm')),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.position='none',
        legend.box = 'horizontal',
        legend.box.just = 'left',
        legend.key.size = unit(0.8, 'cm'))

p5 <-
  ggplot()+
  geom_bar(stat = "identity", 
           data = CN_ANA_LPJ,
           mapping=aes(ModPTRT, deltaCveg.mean/1000, fill=ModPTRT),
           position=position_dodge(preserve = 'single'), 
           alpha = 0.7,
           col="black",
           width = 0.5)+
  geom_errorbar(stat = "identity",
                data = CN_ANA_LPJ,
                mapping=aes(ModPTRT,
                            ymax=deltaCveg.mean/1000+deltaCveg.sd/1000,
                            ymin=deltaCveg.mean/1000-deltaCveg.sd/1000,
                            group=PTRT),
                # position=position_dodge(preserve = 'single', width=0.7),
                linewidth=1,
                width=0.3)+
  scale_fill_manual(name="",
                    values=c("LPJGN"="cyan4",
                             "LPJGP_NOP"="grey70",
                             "LPJGP_MDP"="goldenrod1",
                             "LPJGP_HIP"="seagreen"))+
  scale_x_discrete(name='',
                   labels = c("LPJGN"="CN",
                              "LPJGP_NOP"="aP",
                              "LPJGP_MDP"="MP",
                              "LPJGP_HIP"="HP"))+
  scale_y_continuous(limits = c(-0.05, 0.7),
                     breaks = seq(0, 0.7, 0.2))+
  geom_hline(yintercept=0, lty=1)+
  ylab(expression(Delta * C[veg] * "  ( kg C " * m^-2 * " " * yr^-1 * " )"))+
  ggtitle('')+
  annotate('text', x=0.73, y=0.67, label='(g)', size=10)+
  theme_linedraw() +
  guides(fill = FALSE)+
  theme(text = element_text(family = 'serif'),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.margin = unit(c(0,0.8,0,0.5), 'cm'),
        plot.title = element_text(size=16, face="bold",
                                  hjust = 0.5),
        # axis.text.x=element_text(size=22, vjust = -0.5),
        axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=22),
        axis.title.y=element_text(size=24, vjust = 2.5,
                                  margin = unit(c(0.28,0.28,0.28,0.28), 'cm')),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.position='none',
        legend.box = 'horizontal',
        legend.box.just = 'left',
        legend.key.size = unit(0.8, 'cm'))

p6 <-
  ggplot()+
  geom_bar(stat = "identity", 
           data = CN_ANA_LPJ,
           mapping=aes(ModPTRT, NEP.mean/1000, fill=ModPTRT),
           position=position_dodge(preserve = 'single'), 
           alpha = 0.7,
           col="black",
           width = 0.5)+
  geom_errorbar(stat = "identity",
                data = CN_ANA_LPJ,
                mapping=aes(ModPTRT,
                            ymax=NEP.mean/1000+NEP.sd/1000,
                            ymin=NEP.mean/1000-NEP.sd/1000,
                            group=PTRT),
                # position=position_dodge(preserve = 'single', width=0.7),
                linewidth=1,
                width=0.3)+
  scale_fill_manual(name="",
                    values=c("LPJGN"="cyan4",
                             "LPJGP_NOP"="grey70",
                             "LPJGP_MDP"="goldenrod1",
                             "LPJGP_HIP"="seagreen"))+
  scale_x_discrete(name='',
                   labels = c("LPJGN"="CN",
                              "LPJGP_NOP"="aP",
                              "LPJGP_MDP"="MP",
                              "LPJGP_HIP"="HP"))+
  # scale_y_continuous(limits = c(0, 0.5),
  #                    breaks = seq(0, 0.5, 0.1))+
  geom_hline(yintercept=0, lty=1)+
  ylab(expression(NEP * "  ( kg C " * m^-2 * " " * yr^-1 * " )"))+
  ggtitle('')+
  annotate('text', x=0.73, y=0.61, label='(k)', size=10)+
  theme_linedraw() +
  guides(fill = FALSE)+
  theme(text = element_text(family = 'serif'),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.margin = unit(c(0,0.8,0.5,0.5), 'cm'),
        plot.title = element_text(size=16, face="bold",
                                  hjust = 0.5),
        axis.text.x=element_text(size=28, vjust = -0.5),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=22),
        axis.title.y=element_text(size=24, vjust = 2.5,
                                  margin = unit(c(0.28,0.28,0.28,0.28), 'cm')),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.position='none',
        legend.box = 'horizontal',
        legend.box.just = 'left',
        legend.key.size = unit(0.8, 'cm'))


## comparation of simulated CO2 effects on GPP, biomass growth, and NEP
# ---------------------------------------------------------------------------------------
CN_ANA1_GDAY <- subset(CN_ANA1, ModName == 'A_GDAYP' | ModName == 'I_GDAYN')
CN_ANA1_GDAY$ModPTRT <- factor(CN_ANA1_GDAY$ModPTRT, 
                               levels = c('I_GDAYN_NOP', 'A_GDAYP_NOP', 'A_GDAYP_MDP', 'A_GDAYP_HIP'),
                               labels = c('GDAYN', 'GDAYP_NOP', 'GDAYP_MDP', 'GDAYP_HIP'))
p11 <-
  ggplot()+
  geom_bar(stat = "identity", 
           data = CN_ANA1_GDAY,
           mapping=aes(ModPTRT, GPP.mean/1000, fill=ModPTRT),
           position=position_dodge(preserve = 'single'), 
           alpha = 0.7,
           col="black",
           width = 0.5)+
  geom_errorbar(stat = "identity",
                data = CN_ANA1_GDAY,
                mapping=aes(ModPTRT,
                            ymax=GPP.mean/1000+GPP.sd/1000,
                            ymin=GPP.mean/1000-GPP.sd/1000,
                            group=PTRT),
                # position=position_dodge(preserve = 'single', width=0.7),
                linewidth=1,
                width=0.3)+
  scale_fill_manual(name="",
                    values=c("GDAYN"="cyan4",
                             "GDAYP_NOP"="grey70",
                             "GDAYP_MDP"="goldenrod1",
                             "GDAYP_HIP"="seagreen"))+
  scale_x_discrete(name='',
                   labels = c("GDAYN"="CN",
                              "GDAYP_NOP"="aP",
                              "GDAYP_MDP"="MP",
                              "GDAYP_HIP"="HP"))+
  # scale_y_continuous(limits = c(0, 3),
  #                    breaks = seq(0, 3, 0.5))+
  geom_hline(yintercept=0, lty=1)+
  ylab(expression(CO[2] * " effect on " * GPP * "  ( kg C " * m^-2 * " " * yr^-1 * " )"))+
  ggtitle('GDAY')+
  annotate('text', x=0.73, y=0.78, label='(b)', size=10)+
  theme_linedraw() +
  guides(fill = FALSE)+
  theme(text = element_text(family = 'serif'),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.margin = unit(c(0.5,0.8,0,0.5), 'cm'),
        plot.title = element_text(size=32,
                                  hjust = 0.5),
        # axis.text.x=element_text(size=22, vjust = -0.5),
        axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=22),
        axis.title.y=element_text(size=24, vjust = 2.5,
                                  margin = unit(c(0.28,0.28,0.28,0.28), 'cm')),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.position=c(0.2,0.85),
        legend.box = 'horizontal',
        legend.box.just = 'left',
        legend.key.size = unit(0.8, 'cm'));

p21 <-
  ggplot()+
  geom_bar(stat = "identity", 
           data = CN_ANA1_GDAY,
           mapping=aes(ModPTRT, deltaCveg.mean/1000, fill=ModPTRT),
           position=position_dodge(preserve = 'single'), 
           alpha = 0.7,
           col="black",
           width = 0.5)+
  geom_errorbar(stat = "identity",
                data = CN_ANA1_GDAY,
                mapping=aes(ModPTRT,
                            ymax=deltaCveg.mean/1000+deltaCveg.sd/1000,
                            ymin=deltaCveg.mean/1000-deltaCveg.sd/1000,
                            group=PTRT),
                # position=position_dodge(preserve = 'single', width=0.7),
                linewidth=1,
                width=0.3)+
  scale_fill_manual(name="",
                    values=c("GDAYN"="cyan4",
                             "GDAYP_NOP"="grey70",
                             "GDAYP_MDP"="goldenrod1",
                             "GDAYP_HIP"="seagreen"))+
  scale_x_discrete(name='',
                   labels = c("GDAYN"="CN",
                              "GDAYP_NOP"="aP",
                              "GDAYP_MDP"="MP",
                              "GDAYP_HIP"="HP"))+
  # scale_y_continuous(limits = c(-0.05, 0.4),
  #                    breaks = seq(0, 0.4, 0.1))+
  geom_hline(yintercept=0, lty=1)+
  ylab(expression(CO[2] * " effect on " * Delta * C[veg] * "  ( kg C " * m^-2 * " " * yr^-1 * " )"))+
  ggtitle('')+
  annotate('text', x=0.73, y=0.32, label='(f)', size=10)+
  theme_linedraw() +
  guides(fill = FALSE)+
  theme(text = element_text(family = 'serif'),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.margin = unit(c(0,0.8,0,0.5), 'cm'),
        plot.title = element_text(size=16, face="bold",
                                  hjust = 0.5),
        # axis.text.x=element_text(size=22, vjust = -0.5),
        axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=22),
        axis.title.y=element_text(size=24, vjust = 2.5,
                                  margin = unit(c(0.28,0.28,0.28,0.28), 'cm')),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.position='none',
        legend.box = 'horizontal',
        legend.box.just = 'left',
        legend.key.size = unit(0.8, 'cm')); 

p31 <-
  ggplot()+
  geom_bar(stat = "identity", 
           data = CN_ANA1_GDAY,
           mapping=aes(ModPTRT, NEP.mean/1000, fill=ModPTRT),
           position=position_dodge(preserve = 'single'), 
           alpha = 0.7,
           col="black",
           width = 0.5)+
  geom_errorbar(stat = "identity",
                data = CN_ANA1_GDAY,
                mapping=aes(ModPTRT,
                            ymax=NEP.mean/1000+NEP.sd/1000,
                            ymin=NEP.mean/1000-NEP.sd/1000,
                            group=PTRT),
                # position=position_dodge(preserve = 'single', width=0.7),
                linewidth=1,
                width=0.35)+
  scale_fill_manual(name="",
                    values=c("GDAYN"="cyan4",
                             "GDAYP_NOP"="grey70",
                             "GDAYP_MDP"="goldenrod1",
                             "GDAYP_HIP"="seagreen"))+
  scale_x_discrete(name='',
                   labels = c("GDAYN"="CN",
                              "GDAYP_NOP"="aP",
                              "GDAYP_MDP"="MP",
                              "GDAYP_HIP"="HP"))+
  scale_y_continuous(limits = c(-0.05, 0.5),
                     breaks = seq(0, 0.5, 0.1))+
  geom_hline(yintercept=0, lty=1)+
  ylab(expression(CO[2] * " effect on " * NEP * "  ( kg C " * m^-2 * " " * yr^-1 * " )"))+
  ggtitle('')+
  annotate('text', x=0.73, y=0.47, label='(j)', size=10)+
  theme_linedraw() +
  guides(fill = FALSE)+
  theme(text = element_text(family = 'serif'),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.margin = unit(c(0,0.8,0.5,0.5), 'cm'),
        plot.title = element_text(size=16, face="bold",
                                  hjust = 0.5),
        axis.text.x=element_text(size=28, vjust = -0.5),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=22),
        axis.title.y=element_text(size=24, vjust = 2.5,
                                  margin = unit(c(0.28,0.28,0.28,0.28), 'cm')),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.position='none',
        legend.box = 'horizontal',
        legend.box.just = 'left',
        legend.key.size = unit(0.8, 'cm')); 


## LPJ-GUESS
CN_ANA1_LPJ <- subset(CN_ANA1, ModName == 'D_LPJGP' | ModName == 'J_LPJGN')
CN_ANA1_LPJ$ModPTRT <- factor(CN_ANA1_LPJ$ModPTRT, 
                              levels = c('J_LPJGN_NOP', 'D_LPJGP_NOP', 'D_LPJGP_MDP', 'D_LPJGP_HIP'),
                              labels = c('LPJGN', 'LPJGP_NOP', 'LPJGP_MDP', 'LPJGP_HIP'))
p41 <-
  ggplot()+
  geom_bar(stat = "identity", 
           data = CN_ANA1_LPJ,
           mapping=aes(ModPTRT, GPP.mean/1000, fill=ModPTRT),
           position=position_dodge(preserve = 'single'), 
           alpha = 0.7,
           col="black",
           width = 0.5)+
  geom_errorbar(stat = "identity",
                data = CN_ANA1_LPJ,
                mapping=aes(ModPTRT,
                            ymax=GPP.mean/1000+GPP.sd/1000,
                            ymin=GPP.mean/1000-GPP.sd/1000,
                            group=PTRT),
                # position=position_dodge(preserve = 'single', width=0.7),
                linewidth=1,
                width=0.3)+
  scale_fill_manual(name="",
                    values=c("LPJGN"="cyan4",
                             "LPJGP_NOP"="grey70",
                             "LPJGP_MDP"="goldenrod1",
                             "LPJGP_HIP"="seagreen"))+
  scale_x_discrete(name='',
                   labels = c("LPJGN"="CN",
                              "LPJGP_NOP"="aP",
                              "LPJGP_MDP"="MP",
                              "LPJGP_HIP"="HP"))+
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.2))+
  geom_hline(yintercept=0, lty=1)+
  ylab(expression(CO[2] * ' effect on ' * GPP * "  ( kg C " * m^-2 * " " * yr^-1 * " )"))+
  ggtitle('LPJ-GUESS')+
  annotate('text', x=0.73, y=0.965, label='(d)', size=10)+
  theme_linedraw() +
  guides(fill = FALSE)+
  theme(text = element_text(family = 'serif'),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.margin = unit(c(0.5,0.8,0,0.5), 'cm'),
        plot.title = element_text(size=32,
                                  hjust = 0.5),
        # axis.text.x=element_text(size=22, vjust = -0.5),
        axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=22),
        axis.title.y=element_text(size=24, vjust = 2.5,
                                  margin = unit(c(0.28,0.28,0.28,0.28), 'cm')),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.position='none',
        legend.box = 'horizontal',
        legend.box.just = 'left',
        legend.key.size = unit(0.8, 'cm')); 

p51 <-
  ggplot()+
  geom_bar(stat = "identity", 
           data = CN_ANA1_LPJ,
           mapping=aes(ModPTRT, deltaCveg.mean/1000, fill=ModPTRT),
           position=position_dodge(preserve = 'single'), 
           alpha = 0.7,
           col="black",
           width = 0.5)+
  geom_errorbar(stat = "identity",
                data = CN_ANA1_LPJ,
                mapping=aes(ModPTRT,
                            ymax=deltaCveg.mean/1000+deltaCveg.sd/1000,
                            ymin=deltaCveg.mean/1000-deltaCveg.sd/1000,
                            group=PTRT),
                # position=position_dodge(preserve = 'single', width=0.7),
                linewidth=1,
                width=0.3)+
  scale_fill_manual(name="",
                    values=c("LPJGN"="cyan4",
                             "LPJGP_NOP"="grey70",
                             "LPJGP_MDP"="goldenrod1",
                             "LPJGP_HIP"="seagreen"))+
  scale_x_discrete(name='',
                   labels = c("LPJGN"="CN",
                              "LPJGP_NOP"="aP",
                              "LPJGP_MDP"="MP",
                              "LPJGP_HIP"="HP"))+
  # scale_y_continuous(limits = c(-0.05, 0.7),
  #                    breaks = seq(0, 0.7, 0.2))+
  geom_hline(yintercept=0, lty=1)+
  ylab(expression(CO[2] * ' effect on ' *  Delta  * C[veg] * "  ( kg C " * m^-2 * " " * yr^-1 * " )"))+
  ggtitle('')+
  annotate('text', x=0.73, y=0.62, label='(h)', size=10)+
  theme_linedraw() +
  guides(fill = FALSE)+
  theme(text = element_text(family = 'serif'),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.margin = unit(c(0,0.8,0,0.5), 'cm'),
        plot.title = element_text(size=16, face="bold",
                                  hjust = 0.5),
        # axis.text.x=element_text(size=22, vjust = -0.5),
        axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=22),
        axis.title.y=element_text(size=24, vjust = 2.5,
                                  margin = unit(c(0.28,0.28,0.28,0.28), 'cm')),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.position='none',
        legend.box = 'horizontal',
        legend.box.just = 'left',
        legend.key.size = unit(0.8, 'cm')); 

p61 <-
  ggplot()+
  geom_bar(stat = "identity", 
           data = CN_ANA1_LPJ,
           mapping=aes(ModPTRT, NEP.mean/1000, fill=ModPTRT),
           position=position_dodge(preserve = 'single'), 
           alpha = 0.7,
           col="black",
           width = 0.5)+
  geom_errorbar(stat = "identity",
                data = CN_ANA1_LPJ,
                mapping=aes(ModPTRT,
                            ymax=NEP.mean/1000+NEP.sd/1000,
                            ymin=NEP.mean/1000-NEP.sd/1000,
                            group=PTRT),
                # position=position_dodge(preserve = 'single', width=0.7),
                linewidth=1,
                width=0.3)+
  scale_fill_manual(name="",
                    values=c("LPJGN"="cyan4",
                             "LPJGP_NOP"="grey70",
                             "LPJGP_MDP"="goldenrod1",
                             "LPJGP_HIP"="seagreen"))+
  scale_x_discrete(name='',
                   labels = c("LPJGN"="CN",
                              "LPJGP_NOP"="aP",
                              "LPJGP_MDP"="MP",
                              "LPJGP_HIP"="HP"))+
  scale_y_continuous(limits = c(0, 0.6),
                     breaks = seq(0, 0.6, 0.2))+
  geom_hline(yintercept=0, lty=1)+
  ylab(expression(CO[2] * " effect on " * NEP * "  ( kg C " * m^-2 * " " * yr^-1 * " )"))+
  ggtitle('')+
  annotate('text', x=0.73, y=0.575, label='(l)', size=10)+
  theme_linedraw() +
  guides(fill = FALSE)+
  theme(text = element_text(family = 'serif'),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.margin = unit(c(0,0.8,0.5,0.5), 'cm'),
        plot.title = element_text(size=16, face="bold",
                                  hjust = 0.5),
        axis.text.x=element_text(size=28, vjust = -0.5),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=22),
        axis.title.y=element_text(size=24, vjust = 2.5,
                                  margin = unit(c(0.28,0.28,0.28,0.28), 'cm')),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.position='none',
        legend.box = 'horizontal',
        legend.box.just = 'left',
        legend.key.size = unit(0.8, 'cm'));

plot_comparation <- 
  p1 + p11 + p4 + p41 +
  p2 + p21 + p5 + p51 +
  p3 + p31 + p6 + p61 + 
  plot_layout(ncol = 4)





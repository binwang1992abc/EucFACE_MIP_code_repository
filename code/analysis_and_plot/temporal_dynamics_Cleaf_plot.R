


#############################################################################
#### Modelling CO2 x P at EucFACE: multi-model intercomparison           ####
#### Temporal dynamics of Cleaf with P and CO2 treatments                ####
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


#### data prepare and plot
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

## organize to make neat
ambDF$CO2TRT <- 'ambient'
eleDF$CO2TRT <- 'elevated'

## select 2018-2024
ambDF <- subset(ambDF, YEAR >= 2018 & YEAR <= 2024 & PTRT != 'MDP', select = c('ModName', 'YEAR', 'CO2TRT', 'PTRT', 'CL'))
eleDF <- subset(eleDF, YEAR >= 2018 & YEAR <= 2024 & PTRT != 'MDP', select = c('ModName', 'YEAR', 'CO2TRT', 'PTRT', 'CL'))

## combine 
myDF <- rbind(ambDF, eleDF)

myDF$CP <- factor(paste(myDF$CO2TRT, myDF$PTRT, sep = '_'))
myDF$PTRT <- factor(myDF$PTRT, levels = c('NOP', 'HIP'), labels = c('aP', 'eP'))
myDF$CL[myDF$PTRT == 'eP' & myDF$YEAR < 2020] <- NA

## plot
p1 <-
ggplot(subset(myDF, ModName == 'A_GDAYP'),
       aes(YEAR, CL/1000, group=CP, colour = PTRT))+
  geom_line() +
  geom_point(aes(shape = CO2TRT), size=5) +
  scale_color_manual(names(''),
                     values = c('lightblue3', 'orangered3'))+
  scale_shape_manual(names(''),
                     values = c(22,21),
                     labels = c('ambient' = expression(aCO[2]),
                                'elevated' = expression(eCO[2])))+
  scale_y_continuous(limits = c(0.19,0.38),
                     breaks = seq(0.20,0.375,0.025))+
  annotate('text', x=2018.5, y=0.37, label='(a)', size=10)+
  ggtitle('GDAYP') +
  ylab(expression(C[leaf] * " (kg C " * m^-2 * " " * yr^-1 * ")"))+
  annotate("rect", 
           xmin = 2020, xmax = 2022, 
           ymin = -Inf, ymax = Inf, 
           fill = "grey", alpha = 0.5)+
  theme(
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black', linewidth = 1),
        plot.margin = unit(c(0.5,0.2,0,1.5), 'cm'),
        plot.title = element_text(size=24, 
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
  
p2 <-
  ggplot(subset(myDF, ModName == 'H_QUJSM'),
         aes(YEAR, CL/1000, group=CP, colour = PTRT))+
  geom_line() +
  geom_point(aes(shape = CO2TRT), size=5) +
  scale_color_manual(names(''),
                     values = c('lightblue3', 'orangered3'))+
  scale_shape_manual(names(''),
                     values = c(22,21),
                     labels = c('ambient' = expression(aCO[2]),
                                'elevated' = expression(eCO[2])))+
  scale_y_continuous(limits = c(0.18,0.29),
                     breaks = seq(0.18,0.28,0.02))+
  annotate('text', x=2018.5, y=0.28, label='(b)', size=10)+
  ggtitle('QUJSM') +
  ylab(NULL)+
  annotate("rect", 
           xmin = 2020, xmax = 2022, 
           ymin = -Inf, ymax = Inf, 
           fill = "grey", alpha = 0.5)+
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major=element_blank(),
    panel.background = element_rect(fill = 'white', colour = 'black', linewidth = 1),
    plot.margin = unit(c(0.1,0.2,0,1.5), 'cm'),
    plot.title = element_text(size=24,
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

  p3 <-
  ggplot(subset(myDF, ModName == 'C_CABLP'),
         aes(YEAR, CL/1000, group=CP, colour = PTRT))+
    geom_line() +
    geom_point(aes(shape = CO2TRT), size=5) +
    scale_color_manual(names(''),
                       values = c('lightblue3', 'orangered3'))+
    scale_shape_manual(names(''),
                       values = c(22,21),
                       labels = c('ambient' = expression(aCO[2]),
                                  'elevated' = expression(eCO[2])))+
    scale_y_continuous(limits = c(0.220,0.405),
                       breaks = seq(0.225,0.400,0.025))+
    annotate('text', x=2018.5, y=0.39, label='(c)', size=10)+
    ggtitle('CABLP') +
    ylab(NULL)+
    annotate("rect", 
             xmin = 2020, xmax = 2022, 
             ymin = -Inf, ymax = Inf, 
             fill = "grey", alpha = 0.5)+
    theme(
      panel.grid.minor=element_blank(),
      panel.grid.major=element_blank(),
      panel.background = element_rect(fill = 'white', colour = 'black', linewidth = 1),
      plot.margin = unit(c(0.1,0.2,0,1.5), 'cm'),
      plot.title = element_text(size=24,
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
  
  p4 <-
  ggplot(subset(myDF, ModName == 'E_OCHDP'),
         aes(YEAR, CL/1000, group=CP, colour = PTRT))+
    geom_line() +
    geom_point(aes(shape = CO2TRT), size=5) +
    scale_color_manual(names(''),
                       values = c('lightblue3', 'orangered3'))+
    scale_shape_manual(names(''),
                       values = c(22,21),
                       labels = c('ambient' = expression(aCO[2]),
                                  'elevated' = expression(eCO[2])))+
    scale_y_continuous(limits = c(0.33,0.68),
                       breaks = seq(0.35,0.65,0.05))+
    annotate('text', x=2018.5, y=0.65, label='(d)', size=10)+
    ggtitle('OCHDP') +
    ylab(NULL)+
    annotate("rect", 
             xmin = 2020, xmax = 2022, 
             ymin = -Inf, ymax = Inf, 
             fill = "grey", alpha = 0.5)+
    theme(
      panel.grid.minor=element_blank(),
      panel.grid.major=element_blank(),
      panel.background = element_rect(fill = 'white', colour = 'black', linewidth = 1),
      plot.margin = unit(c(0.1,0.2,0.2,1.5), 'cm'),
      plot.title = element_text(size=24,
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
  
  p5 <-
    ggplot(subset(myDF, ModName == 'G_OCHDX'),
           aes(YEAR, CL/1000, group=CP, colour = PTRT))+
    geom_line() +
    geom_point(aes(shape = CO2TRT), size=5) +
    scale_color_manual(names(''),
                       values = c('lightblue3', 'orangered3'))+
    scale_shape_manual(names(''),
                       values = c(22,21),
                       labels = c('ambient' = expression(aCO[2]),
                                  'elevated' = expression(eCO[2])))+
    scale_y_continuous(limits = c(0.40,0.70),
                       breaks = seq(0.40,0.70,0.05))+
    annotate('text', x=2018.5, y=0.68, label='(e)', size=10)+
    ggtitle('OCHDX') +
    ylab(expression(C[leaf] * " (kg C " * m^-2 * " " * yr^-1 * ")"))+
    xlab('Year')+
    annotate("rect", 
             xmin = 2020, xmax = 2022, 
             ymin = -Inf, ymax = Inf, 
             fill = "grey", alpha = 0.5)+
    theme(
      panel.grid.minor=element_blank(),
      panel.grid.major=element_blank(),
      panel.background = element_rect(fill = 'white', colour = 'black', linewidth = 1),
      plot.margin = unit(c(0.5,0.2,0,1.5), 'cm'),
      plot.title = element_text(size=24, 
                                hjust = 0.5),
      axis.text.x=element_text(size=20),
      # axis.text.x = element_blank(),
      axis.title.x=element_text(size=20),
      axis.text.y=element_text(size=22),
      axis.title.y=element_text(size=22, vjust = 1),
      axis.ticks.length = unit(0.2, 'cm'),
      legend.text=element_text(size=23),
      legend.title=element_text(size=20),
      legend.position='none',
      legend.box = 'horizontal',
      legend.box.just = 'left',
      legend.key.size = unit(1.3, 'cm'))
  
  p6 <-
  ggplot(subset(myDF, ModName == 'B_ELMV1'),
         aes(YEAR, CL/1000, group=CP, colour = PTRT))+
    geom_line() +
    geom_point(aes(shape = CO2TRT), size=5) +
    scale_color_manual(names(''),
                       values = c('lightblue3', 'orangered3'))+
    scale_shape_manual(names(''),
                       values = c(22,21),
                       labels = c('ambient' = expression(aCO[2]),
                                  'elevated' = expression(eCO[2])))+
    scale_y_continuous(limits = c(0.12,0.18),
                       breaks = seq(0.12,0.18,0.01))+
    annotate('text', x=2018.5, y=0.175, label='(f)', size=10)+
    ggtitle('ELMV1') +
    ylab(NULL)+
    xlab('Year')+
    annotate("rect", 
             xmin = 2020, xmax = 2022, 
             ymin = -Inf, ymax = Inf, 
             fill = "grey", alpha = 0.5)+
    theme(
      panel.grid.minor=element_blank(),
      panel.grid.major=element_blank(),
      panel.background = element_rect(fill = 'white', colour = 'black', linewidth = 1),
      plot.margin = unit(c(0.1,0.2,0,1.5), 'cm'),
      plot.title = element_text(size=24,
                                hjust = 0.5),
      axis.text.x=element_text(size=20),
      # axis.text.x = element_blank(),
      axis.title.x=element_text(size=20),
      axis.text.y=element_text(size=22),
      axis.title.y=element_text(size=22, vjust = 1),
      axis.ticks.length = unit(0.2, 'cm'),
      legend.text=element_text(size=23),
      legend.title=element_text(size=20),
      legend.position='none',
      legend.box = 'horizontal',
      legend.box.just = 'left',
      legend.key.size = unit(1.3, 'cm'))
  
  p7 <-
  ggplot(subset(myDF, ModName == 'D_LPJGP'),
         aes(YEAR, CL/1000, group=CP, colour = PTRT))+
    geom_line() +
    geom_point(aes(shape = CO2TRT), size=5) +
    scale_color_manual(names(''),
                       values = c('lightblue3', 'orangered3'))+
    scale_shape_manual(names(''),
                       values = c(22,21),
                       labels = c('ambient' = expression(aCO[2]),
                                  'elevated' = expression(eCO[2])))+
    scale_y_continuous(limits = c(0.12,0.36),
                       breaks = seq(0.15,0.35,0.05))+
    annotate('text', x=2018.5, y=0.34, label='(g)', size=10)+
    ggtitle('LPJGP') +
    ylab(NULL)+
    xlab('Year')+
    annotate("rect", 
             xmin = 2020, xmax = 2022, 
             ymin = -Inf, ymax = Inf, 
             fill = "grey", alpha = 0.5)+
    theme(
      panel.grid.minor=element_blank(),
      panel.grid.major=element_blank(),
      panel.background = element_rect(fill = 'white', colour = 'black', linewidth = 1),
      plot.margin = unit(c(0.1,0.2,0,1.5), 'cm'),
      plot.title = element_text(size=24,
                                hjust = 0.5),
      axis.text.x=element_text(size=20),
      # axis.text.x = element_blank(),
      axis.title.x=element_text(size=20),
      axis.text.y=element_text(size=22),
      axis.title.y=element_text(size=22, vjust = 1),
      axis.ticks.length = unit(0.2, 'cm'),
      legend.text=element_text(size=23),
      legend.title=element_text(size=20),
      legend.position='none',
      legend.box = 'horizontal',
      legend.box.just = 'left',
      legend.key.size = unit(1.3, 'cm'))
  
  p8 <-
  ggplot(subset(myDF, ModName == 'F_QUINC'),
         aes(YEAR, CL/1000, group=CP, colour = PTRT))+
    geom_line() +
    geom_point(aes(shape = CO2TRT), size=5) +
    scale_color_manual(names(''),
                       values = c('lightblue3', 'orangered3'))+
    scale_shape_manual(names(''),
                       values = c(22,21),
                       labels = c('ambient' = expression(aCO[2]),
                                  'elevated' = expression(eCO[2])))+
    scale_y_continuous(limits = c(0.27,0.39),
                       breaks = seq(0.28,0.38,0.02))+
    annotate('text', x=2018.5, y=0.38, label='(h)', size=10)+
    ggtitle('QUINC') +
    ylab(NULL)+
    xlab('Year')+
    annotate("rect", 
             xmin = 2020, xmax = 2022, 
             ymin = -Inf, ymax = Inf, 
             fill = "grey", alpha = 0.5)+
    theme(
      panel.grid.minor=element_blank(),
      panel.grid.major=element_blank(),
      panel.background = element_rect(fill = 'white', colour = 'black', linewidth = 1),
      plot.margin = unit(c(0.1,0.2,0.2,1.5), 'cm'),
      plot.title = element_text(size=24,
                                hjust = 0.5),
      axis.text.x=element_text(size=20),
      # axis.text.x = element_blank(),
      axis.title.x=element_text(size=20),
      axis.text.y=element_text(size=22),
      axis.title.y=element_text(size=22, vjust = 1),
      axis.ticks.length = unit(0.2, 'cm'),
      legend.text=element_text(size=23),
      legend.title=element_text(size=20),
      legend.position='none',
      legend.box = 'horizontal',
      legend.box.just = 'left',
      legend.key.size = unit(1.3, 'cm'))
 
  
  plot_out <-
  (p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8) +
    plot_layout(ncol = 4, byrow = T, guides = 'collect') & 
    theme(legend.position = 'bottom',
          legend.key = element_blank(),          
          legend.box.background = element_blank())
  
  
  
  














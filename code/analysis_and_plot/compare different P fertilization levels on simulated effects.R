

##################################################################################
#### Modelling CO2 x P at EucFACE: multi-model intercomparison                ####
#### Compare P fertilizaion levels on simulated effects， simulated by gday   ####
#### Created by Bin Wang, 2025-04-27                                          ####
#### Email: wangbin1992@zju.edu.cn                                            ####
##################################################################################


rm(list = ls());


## library
library(doBy)
library(ggplot2)
library(cowplot)
library(patchwork)

## color
SpectralPalette <- c('grey70', colorRampPalette(c( 'yellow3', "goldenrod1", 'seagreen'))(10)) 

## set path
setwd("D:/Analysis/EucFACE_MIP/Model_simulation/GDAY_EucFACE_NEW/output_annual_2/")

## read GDAY annual data, new simulation, screening
GDAY_ann <- readRDS(paste0(getwd(), "/data/compile_output/GDAYP_ALL_CO2_x_P_annual.rds"))
GDAY_ann$Cveg <- rowSums(data.frame(GDAY_ann$CL, GDAY_ann$CW, GDAY_ann$CFR, GDAY_ann$CCR, GDAY_ann$CSTOR), na.rm = T)
GDAY_ann$deltaCveg <- rowSums(data.frame(GDAY_ann$deltaCL, GDAY_ann$deltaCW, GDAY_ann$deltaCFR,
                                         GDAY_ann$deltaCCR, GDAY_ann$deltaCSTOR), na.rm = TRUE)

## select P fertilization in the first day of year 2013-2015, fixed climate
myDF <- subset(GDAY_ann, PFERTTYPE == 'day' & CLIM == 'FIX')

## plot P fertilization effect, and its impact on CO2 effect
myDF1 <- subset(myDF, YEAR <= 2015 & YEAR >= 2013)
ambDF <- subset(myDF1, CO2TRT == 'AMB')
ambDF1 <- summaryBy(.~PTRT, data = ambDF[, -c(1,156,157,159)], FUN = c(mean, sd), na.rm = T, keep.names = T)

p1 <-
  ggplot()+
  geom_bar(stat = "identity", 
           data = ambDF1,
           mapping=aes(PTRT, GPP.mean/1000, group=PTRT, fill=PTRT),
           position=position_dodge(preserve = 'single'), 
           alpha = 0.7,
           col="black",
           width = 0.7)+ 
  geom_errorbar(stat = "identity", 
                data = ambDF1,
                mapping=aes(PTRT, 
                            ymax=GPP.mean/1000+GPP.sd/1000,
                            ymin=GPP.mean/1000-GPP.sd/1000,
                            group=PTRT),
                position=position_dodge(preserve = 'single', width=0.7), 
                width=0.35)+
  scale_fill_manual(name="",
                    values=SpectralPalette)+
  scale_y_continuous(limits = c(0, 2.5),
                     breaks = seq(0, 2.5, 0.5))+
  geom_hline(yintercept=0, lty=1)+
  ylab(expression(GPP * " (kg C " * m^-2 * " " * yr^-1 * ")"))+
  ggtitle('') +
  annotate('text', x=1, y=2.4, label='(a)', size=8)+
  theme_linedraw() +
  theme(text = element_text(family = 'serif'),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.margin = unit(c(0.4,0.8,0,0.6), 'cm'),
        plot.title = element_text(size=16, face="bold",
                                  hjust = 0.5),
        # axis.text.x=element_text(size=18, vjust = -0.4),
        axis.text = element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size=20, vjust = 2.5,
                                  margin = unit(c(0,0,0,0), 'cm')),
        legend.text=element_text(size=18),
        legend.title=element_text(size=20),
        legend.position='none',
        legend.box = 'horizontal',
        legend.box.just = 'left',
        legend.key.size = unit(0.8, 'cm'))

p3 <-
  ggplot()+
  geom_bar(stat = "identity", 
           data = ambDF1,
           mapping=aes(PTRT, deltaCveg.mean/1000, group=PTRT, fill=PTRT),
           position=position_dodge(preserve = 'single'), 
           alpha = 0.7,
           col="black",
           width = 0.7)+
  geom_errorbar(stat = "identity", 
                data = ambDF1,
                mapping=aes(PTRT, 
                            ymax=deltaCveg.mean/1000+deltaCveg.sd/1000,
                            ymin=deltaCveg.mean/1000-deltaCveg.sd/1000,
                            group=PTRT),
                position=position_dodge(preserve = 'single', width=0.7), 
                width=0.35)+
  scale_fill_manual(name="",
                    values=SpectralPalette)+
  scale_color_manual(name="",
                     values=SpectralPalette)+
  # scale_y_continuous(limits = c(-0.03, 0.3),
  #                  breaks = seq(0, 0.3, 0.1))+
  geom_hline(yintercept=0, lty=1)+
  ylab(expression(Delta * C[veg] * " (kg C " * m^-2 * " " * yr^-1 * ")"))+
  ggtitle('') +
  guides(fill='none')+
  annotate('text', x=1, y=0.29, label='(c)', size=8)+
  theme_linedraw() +
  theme(text = element_text(family = 'serif'),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.margin = unit(c(0,0.8,0,0.6), 'cm'),
        plot.title = element_text(size=16, face="bold",
                                  hjust = 0.5),
        # axis.text.x=element_text(size=18, vjust = -0.4),
        axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size=20, vjust = 2.5,
                                  margin = unit(c(0,0,0,0), 'cm')),
        legend.text=element_text(size=18),
        legend.title=element_text(size=20),
        legend.position='none',
        legend.box = 'horizontal',
        legend.box.just = 'left',
        legend.key.size = unit(0.8, 'cm'))

p5 <-
  ggplot()+
  geom_bar(stat = "identity", 
           data = ambDF1,
           mapping=aes(PTRT, NEP.mean/1000, group=PTRT, fill=PTRT),
           position=position_dodge(preserve = 'single'), 
           alpha = 0.7,
           col="black",
           width = 0.7)+ 
  geom_errorbar(stat = "identity", 
                data = ambDF1,
                mapping=aes(PTRT, 
                            ymax=NEP.mean/1000+NEP.sd/1000,
                            ymin=NEP.mean/1000-NEP.sd/1000,
                            group=PTRT),
                position=position_dodge(preserve = 'single', width=0.7), 
                width=0.35)+
  scale_fill_manual(name="",
                    values=SpectralPalette)+
  scale_x_discrete(name='',
                   labels = c("NOP"="aP",
                              "P0.1"=expression(P[0.1]),
                              "P0.3"=expression(P[0.3]),
                              "P0.5"=expression(P[0.5]),
                              "P0.7"=expression(P[0.7]),
                              "P0.9"=expression(P[0.9]),
                              "P1.2"=expression(P[1.2]),
                              "P1.5"=expression(P[1.5]),
                              "P2.0"=expression(P[2.0]),
                              "P2.5"=expression(P[2.5]),
                              "P3.0"=expression(P[3.0])))+
  scale_y_continuous(limits = c(-0.02, 0.42),
                     breaks = seq(0, 0.4, 0.1))+
  geom_hline(yintercept=0, lty=1)+
  ylab(expression(NEP * " (kg C " * m^-2 * " " * yr^-1 * ")"))+
  ggtitle('') +
  annotate('text', x=1, y=0.4, label='(e)', size=8)+
  theme_linedraw() +
  theme(text = element_text(family = 'serif'),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.margin = unit(c(0,0.8,0.4,0.6), 'cm'),
        plot.title = element_text(size=16, face="bold",
                                  hjust = 0.5),
        axis.text.x=element_text(size=22, vjust = -0.4),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size=20, vjust = 2.5,
                                  margin = unit(c(0,0,0,0), 'cm')),
        legend.text=element_text(size=18),
        legend.title=element_text(size=20),
        legend.position='none',
        legend.box = 'horizontal',
        legend.box.just = 'left',
        legend.key.size = unit(0.8, 'cm'))



### plot P fertilization on CO2 effect, 2013-2015
# ------------------------------------------------------------------
myDF1 <- subset(myDF, YEAR <= 2015 & YEAR >= 2013)
myDF1 <- myDF1[, c(157,158,1,2:154,160:161)]
ambDF <- subset(myDF1, CO2TRT == 'AMB')
eleDF <- subset(myDF1, CO2TRT == 'ELE')

d1 <- dim(ambDF)[2]
co2dif <- ambDF
co2dif[, c(4:d1)] <- eleDF[, c(4:d1)] - ambDF[, c(4:d1)]
co2pct <- co2dif
co2pct[, c(4:d1)] <- co2dif[, c(4:d1)]/ambDF[, c(4:d1)]*100

co2dif1 <- summaryBy(.~PTRT, data = co2dif[, c(2, 11:158)], FUN = c(mean, sd), na.rm = T, keep.names = T)
co2pct1 <- summaryBy(.~PTRT, data = co2pct[, c(2, 11:158)], FUN = c(mean, sd), na.rm = T, keep.names = T)

p2 <-
  ggplot()+
  geom_bar(stat = "identity", 
           data = co2dif1,
           mapping=aes(PTRT, GPP.mean/1000, group=PTRT, fill=PTRT),
           position=position_dodge(preserve = 'single'), 
           alpha = 0.7,
           col="black",
           width = 0.7)+ 
  geom_errorbar(stat = "identity", 
                data = co2dif1,
                mapping=aes(PTRT, 
                            ymax=GPP.mean/1000+GPP.sd/1000,
                            ymin=GPP.mean/1000-GPP.sd/1000,
                            group=PTRT),
                position=position_dodge(preserve = 'single', width=0.7), 
                width=0.35)+
  scale_fill_manual(name="",
                    values=SpectralPalette)+
  # scale_y_continuous(limits = c(0, 8.3),
  #                    breaks = seq(0, 8, 2))+
  geom_hline(yintercept=0, lty=1)+
  ylab(expression(CO[2] * ' effect on  ' * GPP * "  (kg C " * m^-2 * " " * yr^-1 * ")"))+   
  annotate('text', x=1, y=0.165, label='(b)', size=8)+
  ggtitle('') +
  theme_linedraw() +
  theme(text = element_text(family = 'serif'),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.margin = unit(c(0.4,0.8,0,0.6), 'cm'),
        plot.title = element_text(size=16, face="bold",
                                  hjust = 0.5),
        # axis.text.x=element_text(size=18, vjust = -0.4),
        axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size=20, vjust = 2.5,
                                  margin = unit(c(0,0,0,0), 'cm')),
        legend.text=element_text(size=18),
        legend.title=element_text(size=20),
        legend.position='none',
        legend.box = 'horizontal',
        legend.box.just = 'left',
        legend.key.size = unit(0.8, 'cm'))


p4<-
  ggplot()+
  geom_bar(stat = "identity", 
           data = co2dif1,
           mapping=aes(PTRT, deltaCveg.mean/1000, group=PTRT, fill=PTRT),
           position=position_dodge(preserve = 'single'), 
           alpha = 0.7,
           col="black",
           width = 0.7)+ 
  geom_errorbar(stat = "identity", 
                data = co2dif1,
                mapping=aes(PTRT, 
                            ymax=deltaCveg.mean/1000+deltaCveg.sd/1000,
                            ymin=deltaCveg.mean/1000-deltaCveg.sd/1000,
                            group=PTRT),
                position=position_dodge(preserve = 'single', width=0.7), 
                width=0.35)+
  scale_fill_manual(name="",
                    values=SpectralPalette)+
  scale_y_continuous(limits = c(-0.02, 0.04),
                     breaks = seq(-0.02, 0.04, 0.02))+
  geom_hline(yintercept=0, lty=1)+
  ylab(expression(CO[2] * ' effect on ' * Delta * C[veg] * "  (kg C " * m^-2 * " " * yr^-1 * ")"))+
  ggtitle('') +
  annotate('text', x=1, y=0.038, label='(d)', size=8)+
  theme_linedraw() +
  theme(text = element_text(family = 'serif'),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.margin = unit(c(0,0.8,0,0.6), 'cm'),
        plot.title = element_text(size=16, face="bold",
                                  hjust = 0.5),
        # axis.text.x=element_text(size=18, vjust = -0.4),
        axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size=20, vjust = 2.5,
                                  margin = unit(c(0,0,0,0), 'cm')),
        legend.text=element_text(size=18),
        legend.title=element_text(size=20),
        legend.position='none',
        legend.box = 'horizontal',
        legend.box.just = 'left',
        legend.key.size = unit(0.8, 'cm'))



p6 <-
  ggplot()+
  geom_bar(stat = "identity", 
           data = co2dif1,
           mapping=aes(PTRT, NEP.mean/1000, group=PTRT, fill=PTRT),
           position=position_dodge(preserve = 'single'), 
           alpha = 0.7,
           col="black",
           width = 0.7)+ 
  geom_errorbar(stat = "identity", 
                data = co2dif1,
                mapping=aes(PTRT, 
                            ymax=NEP.mean/1000+NEP.sd/1000,
                            ymin=NEP.mean/1000-NEP.sd/1000,
                            group=PTRT),
                position=position_dodge(preserve = 'single', width=0.7), 
                width=0.35)+
  scale_fill_manual(name="",
                    values=SpectralPalette)+
  scale_x_discrete(name='',
                   labels = c("NOP"="aP",
                              "P0.1"=expression(P[0.1]),
                              "P0.3"=expression(P[0.3]),
                              "P0.5"=expression(P[0.5]),
                              "P0.7"=expression(P[0.7]),
                              "P0.9"=expression(P[0.9]),
                              "P1.2"=expression(P[1.2]),
                              "P1.5"=expression(P[1.5]),
                              "P2.0"=expression(P[2.0]),
                              "P2.5"=expression(P[2.5]),
                              "P3.0"=expression(P[3.0])))+
  scale_y_continuous(limits = c(-0.02, 0.067),
                     breaks = seq(-0.02, 0.06, 0.02))+
  geom_hline(yintercept=0, lty=1)+
  ylab(expression(CO[2] * ' effect on NEP ' * " (kg C " * m^-2 * " " * yr^-1 * ")"))+
  ggtitle('') +
  annotate('text', x=1, y=0.064, label='(f)', size=8)+
  theme_linedraw() +
  theme(text = element_text(family = 'serif'),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.margin = unit(c(0,0.8,0.4,0.6), 'cm'),
        plot.title = element_text(size=16, face="bold",
                                  hjust = 0.5),
        axis.text.x=element_text(size=22, vjust = -0.4),
        axis.title.x=element_blank(),
        axis.text.y=element_text(size=20),
        axis.title.y=element_text(size=20, vjust = 2.5,
                                  margin = unit(c(0,0,0,0), 'cm')),
        legend.text=element_text(size=18),
        legend.title=element_text(size=20),
        legend.position='none',
        legend.box = 'horizontal',
        legend.box.just = 'left',
        legend.key.size = unit(0.8, 'cm'))

plot_comparation <-
  p1 + p2 +
  p3 + p4 +
  p5 + p6 +
  plot_layout(widths = c(2,2,2,2,2,2),
              ncol = 2)


# qualitative plots ----
library(here)
library(tidyverse)
library(wesanderson)
library(ggpubr)

## read attribute data ----
data <- read.csv("perdiz.csv", header = TRUE, as.is = TRUE)
 
raw <- data$raw.mat # raw material
con <- data$context # burial context
temp <- data$temporal # temporal period
site <- data$trinomial # site

# barplots qual ----
# barplot of raw material count by burial context
raw.con <- ggplot(data, aes(raw)) +
  geom_bar(aes(fill = con))+
  scale_fill_manual(values = wes_palette("Moonrise2")) +
  labs(x = 'Raw Material', y = 'Count')

# barplot of raw material count by temporal
raw.temp <- ggplot(data, aes(raw)) +
  geom_bar(aes(fill = temp))+
  scale_fill_manual(values = wes_palette("Moonrise2")) +
  labs(x = 'Raw Material', y = 'Count')

# barplot of raw material count by site
raw.site <- ggplot(data, aes(trinomial)) +
  geom_bar(aes(fill = raw))+
  scale_fill_manual(values = wes_palette("Moonrise2")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = 'Site', y = 'Count')

# render figure
raw.figure <- ggarrange(raw.site,
                        ggarrange(raw.con,raw.temp,
                                  ncol = 2, 
                                  labels = c("b","c")),
                        nrow = 2,
                        labels = c("a")
)

# plot figure
raw.figure

# temporal plots ----
# load ggplot2
library(ggplot2)
library(wesanderson)

# temporal attributes
pal <- wes_palette("Moonrise2", 20, type = "continuous")

# gantt chart of relative dates for perdiz arrow points
temp<-data.frame(Site = c('41SM442','41AN51, Pace McDonald','41AN115',
                          '41CP5, Tuck Carpenter','41CP12, Johns','41CP20, 
                          BJ Horton','41CP220, Kitchen Branch','41CP495, 
                          Sam D. Carpenter Bottom','41HS15, Pine Tree Mound',
                          '41HS235, Keasler','41HS269, C. D. Marsh',
                          'Hickory Creek #2','41NA49, Washington Square Mound', 
                          '41NA206, Spradley', '41SA135, Jack Walton','41SM55, 
                          Bryan Hardy','41SM193, Redwine','41SM195, Wolf',
                          '41SY43, Old Timers','41SY280, Sybs Site'),
                 Date_Range_CE = c(900,1250,900,1430,1430,1500,1400,1500,
                                   1300,1500,1250,1400,1238,1680,1200,1297,1300,
                                   1315,1400,1400), # in years CE
                 end = c(1200,1450,1680,1500,1600,1550,1499,1550,1650,1680,1450,
                         1499,1445,1725,1400,1391,1454,1440,1680,1590) # in years CE
)

# reorder types by beginning of relative date range
temp$Site <- factor(temp$Site, levels = temp$Site[order(temp$Date_Range_CE)])
# arrange figure
type.time <- ggplot(temp, 
                    aes(x = Date_Range_CE, 
                        xend = end, 
                        y = factor(Site,
                                   levels = rev(levels(factor(Site)))), 
                        yend = Site, 
                        color = Site)) +
  geom_segment(size = 2.5) +
  scale_colour_manual(values = pal) +
  theme(legend.position = "none") +
  labs(y = "Site", x = "Date Range CE")

# render figure
type.time

# elliptical fourier analysis ----
# load packages
devtools::install_github("MomX/Momocs")
library(here)
library(Momocs)

# read images
jpg.list <- list.files(here("img.perdiz"), full.names = TRUE)

# read attribute data
att.data <- read.csv("perdiz.csv", header = TRUE, as.is = TRUE)

# attribute to factor
att.data$region <- as.factor(att.data$region)
att.data$temp.reg <- as.factor(att.data$temp.reg)

# generate outlines
outlines <- jpg.list %>%
  import_jpg()

# add attributes
data.out <- Out(outlines, 
                fac = att.data)

# scale, align, rotate, and center specimens
norm.outlines <- data.out %>% 
  coo_scale() %>%
  coo_align() %>%
  coo_rotate() %>% 
  coo_center()

pile(norm.outlines)

# calibrate how many harmonics needed
calibrate_harmonicpower_efourier(norm.outlines, 
                                 nb.h = 30)

# 11 harmonics needed to capture 99.9 percent of variation
calibrate_reconstructions_efourier(norm.outlines, 
                                   range = 1:11)

# generate efa outlines with 11 harmonics
efa.outlines <- efourier(norm.outlines, 
                         nb.h = 11, 
                         norm = TRUE)

# use efa.outlines for pca
pca.outlines <- PCA(efa.outlines)


## efa exploratory analysis ----
# pca 
scree_plot(pca.outlines)

# plot pca by region
plot_PCA(pca.outlines, 
         morphospace_position = "range_axes",
         palette = pal_qual_solarized,
         chullfilled = TRUE,
         ~region,
         axesnames = TRUE,
         morphospace = TRUE,
         eigen = TRUE,
         center_origin = TRUE,
         zoom = 1.25)

# plot pca by temporal + region
plot_PCA(pca.outlines, 
         morphospace_position = "range_axes",
         palette = pal_qual_solarized,
         chullfilled = TRUE,
         ~temp.reg,
         axesnames = TRUE,
         morphospace = TRUE,
         eigen = TRUE,
         center_origin = TRUE,
         zoom = 1.4)


# mean shape + 2sd for the first 10 pcs
PCcontrib(pca.outlines, nax = 1:5)

## efa confirmatory analysis ----
# manova

# shape difference between temporal periods?
MANOVA(pca.outlines, 'region')

# shape difference between temp.reg?
MANOVA(pca.outlines, 'temp.reg')
# which differ?
MANOVA_PW(pca.outlines, 'temp.reg')

## mean shapes ----

# mean shapes

# region
ms.1 <- MSHAPES(efa.outlines, ~region)
plot_MSHAPES(ms.1, size = 0.75)

# temporal + region
ms.1 <- MSHAPES(efa.outlines, ~temp.reg)
plot_MSHAPES(ms.1, size = 0.75)

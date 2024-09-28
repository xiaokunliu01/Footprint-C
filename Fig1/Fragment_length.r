setwd("/Footprint-C/fragment_Length/")
library(RColorBrewer)
library(ggplot2)
library(ggpubr)

png("fragment_Length.png", units = "in", width = 2.5, height = 2.6, res = 300, bg = "transparent")
# tiff("fragment_Length_WT.tiff", units = "in", width = 2.5, height = 2.6, res = 300, bg = "transparent")
# pdf("fragment_Length_dTAG2h.pdf", width = 2.5, height = 2.6, bg = "transparent")


bin_num <- 16

df <- read.table("Footprint-C.length")


p <- 
  ggplot(df,aes(df[,1],df[,2])) +
  geom_bin2d(aes(fill = after_stat(density)), bins = bin_num,na.rm=TRUE) + 
  scale_fill_gradient2(low = "white",high = "black",na.value = "black",
                       # limits=c(0,0.2),
                        breaks=c(0.01, 0.05)
                       ) +
  xlab('Up fragment (bp)') + 
  ylab('Down fragment (bp)') +
  labs(title = "") +
  scale_x_continuous(limits=c(0,160), breaks=seq(0,150,30),expand = expansion(mult = c(0, .03))) + 
  scale_y_continuous(limits=c(0,160), breaks=seq(0,150,30),expand = c(0,0)) +
  guides(fill = guide_colorbar(
    title = "",
    # title.position = "top",
    # title.vjust = 0.1,
    # label.position = "top",
    label.vjust = 3,
    # label.hjust = -1,
    barwidth = 4,
    barheight = 0.5,
    # nbin=10,
    # raster=TRUE,
    direction="horizontal",
    label.theme = element_text(
                               size = 10,
                               colour = "black",
                               # angle = 270
                               ),
    # keywidth = 0,
    reverse = FALSE
  )) +
  theme(axis.text = element_text(color = "black", size = 10),#face = "bold",
        axis.title = element_text(color = "black", size = 10),#face = "bold",
        panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        # legend.title = element_blank(),
        legend.position = c(0.66,1.05),
        legend.text = element_text(size = 10),#face = "bold"
        legend.title = element_text(color = "black",size = 10),
        rect = element_rect(fill = NA, size=0, color = NA),
        plot.background = element_rect(fill = NA, size=0, color = NA),
        legend.background = element_rect(fill = NA, size=0, color = NA),
        legend.key = element_rect(fill = NA, size=0, color = NA),
  )

p

dev.off()

# plot heatmap based motifs

setwd(dir = "/home/lxk/private/Footprint-C/FootprintC_paper/heatmap_CTCF_CTCF/3w/CTCF_MAZ_CTCF/121275_R1/")

# png("tmp4.png", units = "in", width = 10, height = 10, res = 300, bg = "transparent")
pdf("RR_2.pdf", width = 3.7, height = 3.7, bg = "transparent")

df <- read.table("RR_2.txt",sep = "\t", header = F)
p <- df %>%
  ggplot(aes(x = as.character(df[,1]), y = as.character(df[,2]))) +
  geom_tile(aes(fill = df[,3])) +
  scale_fill_gradient2(low = "white",
                       # high = "firebrick",
                       #B22222, FR
                       # high = "springgreen4",
                       #008B45, FF
                       # high = "blue3",
                       #0000CD, RF
                       high = "turquoise3",
                       #00C5CD, RR
                       # firebrick springgreen4 blue3 turquoise3
                       na.value = "#DCDCDC",
                       limit =c(0,22),
                       breaks=c(0,7,14,21)
  ) +
  guides(fill = guide_colorbar(
    title = "",
    label.hjust = unit(-0.1,"cm"),
    # label.theme = element_text(
    #   size = 10,
    #   colour = "black",
    # )
  )
  ) +
  xlab("") + # up anchor
  ylab("") + # down anchor
  # labs(title = "CTCF HiChIP Tandem F loop") +
  scale_x_discrete(
    limits = c("121271","121272","121273","121274","121275","121276","121277","121278"),
    breaks=c("121271","121272","121273","121274","121275","121276","121277","121278"),
    labels=c("M53183","C68089","M53184","M53185","C68090","M53186","M53187","C68091"),
    expand = c(0,0)
  ) +
  scale_y_discrete(
    limits = c("121278","121277","121276","121275","121274","121273","121272","121271"),
    breaks=c("121278","121277","121276","121275","121274","121273","121272","121271"),
    labels=c("C68091","M53187","M53186","C68090","M53185","M53184","C68089","M53183"),
    expand = c(0,0)
  ) +
  theme(axis.text = element_text(color = "black", size = 9), 
        # face = "bold",
        axis.title = element_text(color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        # legend.position = "none",
        legend.key.height = unit(0.4, "cm"), 
        legend.key.width = unit(0.3, "cm"),
        legend.text = element_text(color = "black",size = 9),
        legend.title = element_text(color = "black",size = 9),
        rect = element_rect(fill = NA, size=0, color = NA),
        plot.background = element_rect(fill = NA, size=0, color = NA),
        legend.background = element_rect(fill = NA, size=0, color = NA),
        legend.key = element_rect(fill = NA, size=0, color = NA),
  ) +
  coord_equal(ratio = 1)



p
dev.off()

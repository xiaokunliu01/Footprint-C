setwd(dir = "/home/lxk/ad1/DNaseC/DNaseC_paper/response/ZNF143/dimerization/K562/distance_density/homo-homo/")

# png("heatmap_cis_dis1k_17list.png", units = "in", width = 10, height = 10, res = 300, bg = "transparent")
pdf("homo_0.1_-KLF.pdf", width = 4, height = 2.5, bg = "transparent")

df <- read.table("homo.tab",sep = "\t", header = F)
df <- df %>% filter(df[,2]!="KLF-KLF")
p <- df %>% 
  ggplot(aes(x = log10(df[,1]))) +
  geom_line( stat = "density", aes(color=df[,2]), size=1.2,)+
  # geom_density(aes(color=df[,2]), size=1.0) +
  scale_color_manual(
    limits = c("ZNF143-ZNF143","EGR1-EGR1","CTCF-CTCF","MAZ-MAZ"),
    values = c("ZNF143-ZNF143"=rgb(0.5,0.78,0.5),
               "EGR1-EGR1"=rgb(0.62,0.47,0.83),
               "CTCF-CTCF"="#386CAF",
               "MAZ-MAZ"=rgb(1,0.6,0.2))
  ) +
  # "#9f79d3"
  guides(color = guide_legend(
    title = "",
    label.hjust = unit(0,"cm"),
    keywidth = unit(0.3,"cm"),
    # keyheight = unit(0, "cm"),
    ncol = 1,
    label.theme = element_text(size = 13)
  ))+
  xlab("Distance (kb)") + 
  ylab("Density") + 
  scale_x_continuous(
    limits = c(2, 6),
    breaks = c(2, 3, 4, 5, 6),
    labels = c(0.1, 1, 10, 100, 1000),
    expand = expansion(mult = c(0, .05))
  ) +
  scale_y_continuous(
    limits = c(0,0.91),
    breaks = c(0.3, 0.6, 0.9),
    labels = c(0.3, 0.6, 0.9),
    expand = c(0,0)
  ) +
  theme(axis.text = element_text(color = "black", size = 15), 
        # face = "bold",
        axis.title = element_text(color = "black", size = 15),
        # axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_blank(), 
        legend.position = c(0.53,0.82),
        # legend.key.width = unit(0, "cm"),
        legend.text = element_text(color = "black",size = 15),
        legend.title = element_text(color = "black",size = 15),
        rect = element_rect(fill = NA, size=0, color = NA),
        plot.background = element_rect(fill = NA, size=0, color = NA),
        legend.background = element_rect(fill = NA, size=0, color = NA),
        legend.key = element_rect(fill = NA, size=0, color = NA),
  )
  
  p
  
dev.off()
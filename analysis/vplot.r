# plot V-plot

setwd(dir = "Footprint-C/vplot/")

library(tidyverse)
library(colorspace)
png("FragmentLengthVsDistance.png", units = "in", width = 2.8, height = 2.1, res = 300, bg = "transparent")
# pdf("FragmentLengthVsDistance.pdf", width = 5, height = 3, bg = "transparent")

df <- read.table("FragmentLengthVsDistance.csv", sep = ";", header = F)
# df_line1 <- data.frame(x = c(0, -40), y = c(0, 75))
# df_line2 <- data.frame(x = c(0, 40), y = c(0, 75))
# df_line3 <- data.frame(x = c(-53.33, 53.33), y = c(0, 0))

p <- ggplot(df, aes(x=df[,2], y=df[,1])) + 
  stat_density_2d(geom = "tile", aes(fill = ..density..^1), contour = FALSE, n=200) +
  scale_fill_gradient2(low = "white", high = "dodgerblue4") +
  # geom_point(alpha = 1/20, size = 1/1, color = "dodgerblue4") +
  # geom_abline(intercept = 0, slope = 1.875, color="red", size = 1.2) +
  # geom_abline(intercept = 0, slope = -1.875, color="red", size = 1.2) +
  # geom_vline(xintercept = -53.33, color="red", size = 1.2) + 
  # geom_vline(xintercept = 53.33, color="red", size = 1.2) +
  # annotate("text", x=0, y=95, label="43.41%", colour = "red") +
  xlab("Distance from CTCF motif (bp)") + ylab("Fragment Length (bp)") + 
  labs(title = "") +
  scale_x_continuous( 
                      limits = c(-180,180),
                      breaks=seq(-180,180,60),
                      expand = expansion(mult = c(0, .03))
                      ) +
  scale_y_continuous(
                     limits = c(0, 160),
                     breaks=seq(0,150,50),
                     expand = c(0,0)
                     ) +
  theme(
        panel.border = element_blank(), 
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(color = "black", size = 10), #face = "bold",
        axis.title = element_text(color = "black", size = 10),
        legend.position = "none",
        legend.title = element_blank(), 
        legend.text = element_text(size = 10),
        rect = element_rect(fill = NA, size=0, color = NA),
        plot.background = element_rect(fill = NA, size=0, color = NA),
        legend.background = element_rect(fill = NA, size=0, color = NA),
        legend.key = element_rect(fill = NA, size=0, color = NA),
        )
p
dev.off()



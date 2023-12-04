
library(reshape2, quietly = T)
library(dplyr, quietly = T)
library(ggplot2, quietly = T)
library(ggsci, quietly = T)
library(cowplot, quietly = T)
library(pheatmap, quietly = T)
library(ggplotify, quietly = T)
library(RColorBrewer, quietly = T)
#library("export")
library(ggsci, quietly = T)
library(data.table)

my_theme2 <- theme(text = element_text(color = "black",  size=14),
        axis.text = element_text(color = "black",  size=12),
        axis.line = element_line(size = 0.6),
        axis.ticks.length  = unit(1, 'mm'),
        axis.ticks = element_line(size = 0.6),
        axis.text.x = element_text(vjust=0.5),
        plot.title = element_text(size=12))

heat <- function(x,y,z){
        m <- x
        m[m==0] <- ""
        pheatmap(x,  
         cluster_rows = F, 
         cluster_cols = F,
         display_numbers = as.matrix(m),
         fontsize_number = 12,
         fontsize = 12,
         fontfamily= "sans",
         number_format = "%.0f",
         main = y,
         legend = F,
         color = c(colorRampPalette(colors = c("white","brown"))(30),colorRampPalette(colors = c("brown","red4"))(10)),
         legend_breaks=z
                )
        } 


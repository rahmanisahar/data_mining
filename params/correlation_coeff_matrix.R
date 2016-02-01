pdf('~/Desktop/coor_plot_for_2by2_1_6_2_col.pdf')
m31 <- read.csv('~/Desktop/project/data_mining/m31/ascii_tables/m31_table_without_UBVRIJHKs.csv',heade=TRUE, sep=',')
small_m31=m31[c(1,2,6),]
small_m31_1=small_m31[c(4:46)]
small_m31_2=small_m31_1[-1]
small_m31_2=small_m31_2[-4]
x_cor=cor(small_m31_2, use='pairwise',method="pearson")
library(ggplot2)
library(corrplot)
corrplot(x_cor, method = "color", insig = "pch", addrect = 3,tl.cex = 0.8,tl.col='black')
dev.off()
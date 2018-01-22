
###  Immunohistochemistry heatmaps ####

library(lattice); library(gplots)

Z <- as.matrix(read.csv2("/Volumes/LETICIA/0_Autoproject/TMA_at_al/DECEMBER_2017/working_IHC_plot_Dic2017.csv", header=T, row.names=1))
head(Z)
y <- Z[,c(2:6)] ## markers with absolute counts
s <- Z[,c(9:16)] # %%  markers with % values
L <- Z[, c(7:8)]  # PDL1 tumor and lymphocytes

pdf(file="..../IHC.pdf", width=30, height=30)
x1 <- levelplot(y, col.regions=colorpanel(40,  "grey", "yellow",  "darkblue"), main="", cexRow=0.3, cexCol=0.3, scales=list(y=list(rot=45), x=list(rot=45)))
x2 <- levelplot(s, col.regions=colorpanel(40,"grey", "yellow",  "darkblue"), main="", cexRow=0.3, cexCol=0.3, scales=list(y=list(rot=45), x=list(rot=45)))
x3 <- levelplot(L, col.regions=colorpanel(40,"grey", "lightgreen",  "purple"), main="", scales=list(y=list(rot=45), x=list(rot=45)))

library(gridExtra)
grid.arrange(x1, x2, x3, ncol=1)
dev.off()

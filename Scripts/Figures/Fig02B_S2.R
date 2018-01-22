###  Oncoprint ###

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
  },
  white = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "white", col = NA))
  },
  MISSENSE = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h*0.66, gp = gpar(fill = "#008000", col = NA))
  },
  TRUNC = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h*0.33, gp = gpar(fill = "black", col = NA))
  },
  INFRAME = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h-unit(0.4, "mm"), gp = gpar(fill = "brown", col = NA))
  },
  OTHER = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h*0.2, gp = gpar(fill = "purple", col = NA))
  },
  GERM = function(x, y, w, h) {
    grid.rect(x, y, w-unit(1, "mm"), h-unit(0.4, "mm"), gp = gpar(fill = "darkblue", col = NA))
  }
)
col = c("MISSENSE" = "#008000", "TRUNC" = "black", "INFRAME" = "brown", "OTHER" = "purple","GERM" = "darkblue", "white"= "white")

## For breast cancer drivers and non-drivers
X <- read.delim("/BreastcancerDrivers.txt", row.names = 1, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = T) 
X = t(as.matrix(X))
library(ComplexHeatmap)
oncoPrint(X, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, 
          column_title = " Non Breast Cancer Drivers",
          show_column_names = TRUE,
          remove_empty_columns = T,
          row_order = NULL, column_order = NULL ,
          heatmap_legend_param = list(title = "Alterations", at = c( "MISSENSE", "TRUNC", "INFRAME", "OTHER" ,"GERM"), 
                                      labels = c("Missense", "Truncating", "In-frame", "Other", "Germline")))

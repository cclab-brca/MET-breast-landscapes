## Heterogeneity barplots ###

ljsp <-read.table(".../all-filtered.maf.txt", sep="\t",header=TRUE, fill=TRUE,comment.char="#",quote="",blank.lines.skip=T, row.names = NULL, stringsAsFactors = FALSE)
ljsp$status <- factor(ljsp$status, levels = c( "private", 'clade', 'stem'), ordered = TRUE)
cbPalette <- c( "darkgreen", "tomato",  "darkblue")
ljsp$mets <- as.character(ljsp$mets)
ljsp$mets <- factor(ljsp$mets, levels=unique(ljsp$mets))
p <- ggplot(ljsp, aes(mets, fill=status)) + geom_bar() 
p + theme(axis.text.x = element_text(face="bold", size=14, angle = 90, hjust = 1), axis.text.y = element_text(size=14)) +  labs(y ="Number of Somatic Mutations") + scale_fill_manual(values=cbPalette) 


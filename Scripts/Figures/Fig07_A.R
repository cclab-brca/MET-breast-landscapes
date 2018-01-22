###### Figure 7A: TCR proportion of reads ### 

library(reshape)
library(ggplot2)
library(ggsignif)
library(cowplot)

### Alpha proportion of reads ###
alpha <- read.delim("Per_read_public_removed_let.txt")
alpha <- alpha[alpha$chain== "alpha",]
alpha_melt <- melt(alpha, id=c("FINAL_ID_Metastasis_Code","Patient"))
alpha_melt$Percent <- alpha_melt$value/100

pdf(file="Per_proportion_reads_seq_public_removed_alpha_variableV3_UPDATED.pdf", width=10, height=8)
my_pallete <- c("darkblue","tomato", "darkgreen")
alpha_melt$variable <- factor(alpha_melt$variable, levels = c( "stem_CDR3", 'clade_CDR3', 'private_CDR3'), ordered = TRUE)
a <- ggplot(alpha_melt, aes(x=as.factor(variable), y=Percent, fill=factor(variable))) +  geom_boxplot()+  geom_signif(comparisons = list(c("private_CDR3", 'stem_CDR3')), 
                                                                                                                      map_signif_level=F, test = "t.test")
a +  scale_y_continuous(breaks = seq(0, 1.00, by = 0.2), limits=c(0,1)) + background_grid(major = 'xy', minor = "xy") + ggtitle("T cell receptor proportion of reads: alpha chain")  + scale_fill_manual(values=my_pallete, name="status") +  labs(x ="" ) + labs(y ="TCR proportion of reads (%)") + theme(axis.text.x = element_text(face="bold", size=14, angle = 90, hjust = 1), axis.text.y = element_text(size=14))
dev.off()


### Beta proportion of reads ###
beta <- read.delim("Per_read_public_removed_let.txt")
beta <- beta[beta$chain== "beta",]
beta <- beta[, -c(3, 7:10)]
head(beta)
unique(beta$FINAL_ID_Metastasis_Code)

beta_melt <- melt(beta, id=c("FINAL_ID_Metastasis_Code","Patient"))
head(beta_melt)
beta_melt$Percent <- beta_melt$value/100

beta_melt <- read.delim("TCR_proportion_reads_beta_chain_UPDATED.txt")

my_pallete <- c("darkblue","tomato", "darkgreen")
beta_melt$variable <- factor(beta_melt$variable, levels = c("stem_CDR3", 'clade_CDR3', 'private_CDR3'), ordered = TRUE)
pdf(file="Per_proportion_reads_seq_public_removed_beta_variable_UPDATED.pdf", width=10, height=8)
a <- ggplot(beta_melt, aes(x=as.factor(variable), y=Percent, fill=factor(variable))) +  geom_boxplot() + geom_signif(comparisons = list(c("clade_CDR3", 'stem_CDR3')), map_signif_level=F, test = "t.test")
a +  scale_y_continuous(breaks = seq(0, 1.00, by = 0.2), limits=c(0,1)) + background_grid(major = 'xy', minor = "xy") + ggtitle("T cell receptor proportion of reads: beta chain")  + scale_fill_manual(values=my_pallete, name="status") +  labs(x ="" ) + labs(y ="TCR proportion of reads (%)") + theme(axis.text.x = element_text(face="bold", size=14, angle = 90, hjust = 1), axis.text.y = element_text(size=14))
dev.off()







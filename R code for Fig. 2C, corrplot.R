library(psych)
library(corrplot)
library(ggplot2)

###Load data
mydata <- read.csv('Table_S9.csv', head=T, skip=1)  #Feature table
cli40 <- read.csv('Table_S1.csv', row.names=1,  skip=1, head=T)  #Metadata (Clinical data)

###Pre-progress the data
mydata[1:5, 1:10]
genera <- mydata[,-c(1,2,3,4,6)]  #Delete taxonomy except Genus
genera_sum <- aggregate(genera[ ,-1], by=list(genera[ ,1]), sum)
rownames(genera_sum) <- genera_sum[, 1]
otu_g <- genera_sum[,-1]
g <- data.frame(t(otu_g))
g_per <- g/rowSums(g)
rowSums(g_per)

###Grouping factor
group_letter <- substring(rownames(g_per), 1, 1)
g_per$group <- ifelse(group_letter=="P", "JIA", "Control")
g_per$group <- as.factor(g_per$group)
str(g_per)

###Select the four genera which differed between the two groups
g_jia <- subset(g_per, group=="JIA")
g_jia_sig <- g_jia[,c("Anaerostipes", "Dialister", "Lachnospira", "Roseburia")]

###The clinical data is merged with the feature table
cli40_jia <- subset(cli40, Groups=='JIA')
cli39 <- cli40_jia[!(rownames(cli40_jia)=='PBJ70'), !(colnames(cli40_jia)==c('Groups', 'Sex', 'subtype'))]
g_cli <- merge.data.frame(g_jia_sig, cli39, by="row.names")
rownames(g_cli) <- g_cli[, 1]
d <- g_cli[,-1]

###Correlation analysis
aa <- corr.test(d, method="spearman", adjust="holm")
aa_r <- data.frame(aa$r)
aa_p <- data.frame(aa$p)
#Rows with P-values ＜0.05
name_A <- rownames(aa_p)[aa_p$Anaerostipes < 0.05]
name_D <- rownames(aa_p)[aa_p$Dialister < 0.05]
name_L <- rownames(aa_p)[aa_p$Lachnospira < 0.05]
name_R <- rownames(aa_p)[aa_p$Roseburia < 0.05]

##Visualization
cli_1 <- c(name_A, name_D, name_L, name_R)
cli_2 <- unique(cli_1)
cli_3 <- cli_2[5:15]
plot_r <- aa$r[cli_3, 1:4]
str(plot_r)
corrplot(plot_r, method="pie")
plot_p <- aa$p[cli_3, 1:4]
corrplot(plot_r, method="pie", p.mat=plot_p, sig.level=.05, insig="label_sig", cl.ratio=0.5, cl.cex=1.3, 
         tl.cex=1.3, tl.col="black", tl.srt=30)

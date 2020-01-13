library(ALDEx2)
library(ggplot2)
library(ggrepel)

#Load and pre-process the data. 
mydata <- read.csv('Table_S9.csv', skip=1, header=T)
mydata[1:5, 1:10]
myotu <- mydata[,-c(1:5)]
rownames(myotu) <- myotu[, 1]
otu <- myotu[,-1]

#Grouping factor
group_letter_otu <- substring(colnames(otu), 1, 1)
conds_otu <- ifelse(group_letter_otu == "P", "JIA", "Control")

#Centred Log-Ratio transformation and Welch’s t test
set.seed(1000)
otu.log <- aldex.clr(otu, conds_otu, mc.samples=128, verbose=TRUE)
otu.test <- aldex.ttest(otu.log, conds_otu, paired.test=FALSE)
head(otu.test)

#Estimate Effect Size
otu.effect <- aldex.effect(otu.log, conds_otu, include.sample.summary=FALSE, verbose=FALSE)

#Merge all data into one object
otu.p.all <- data.frame(otu.test, otu.effect)
head(otu.p.all)

#Sellcet the variables with P values ＜0.05
subset(otu.p.all, otu.p.all$we.ep < 0.05)
subset(otu.p.all, otu.p.all$we.eBH < 0.05)
subset(otu.p.all, otu.p.all$wi.ep < 0.05)
subset(otu.p.all, otu.p.all$wi.eBH < 0.05)
subset(otu.p.all, otu.p.all$effect < -0.5)

otu.t <- data.frame(t(otu))
otu.pern <- otu.t/rowSums(otu.t)
rowSums(otu.pern)
otu.pern$group <- conds_otu
otu.control <- subset(otu.pern, group=='Control')
otu.jia <- subset(otu.pern, group=='JIA')
(apply(otu.jia[,-1391], 2, median)*100)[c("X361727", "X581003", "X369429", "X470382", "X368261")]
(apply(otu.control[,-1391], 2, median)*100)[c("X361727", "X581003", "X369429", "X470382", "X368261")]

otu.p.all$p.effect <- ifelse(otu.p.all$wi.ep < 0.05 & otu.p.all$effect < 0, "down", 
                             ifelse(otu.p.all$wi.ep < 0.05 & otu.p.all$effect > 0, "up", "unchangeed"))
otu.p.all$sig.otu <- ifelse(otu.p.all$effect < -0.5, rownames(otu.p.all)[otu.p.all$effect < -0.5], "")

#Bland-Altman Plot
aldex.plot(otu.p.all, type="MA", test="welch", cutoff=0.15, all.cex=0.7, called.cex=1.1,
           rare.col="grey", called.col="red")
#Volcano
plot(otu.p.all$effect, otu.p.all$wi.ep, log="y", pch=19, main="Effect",
     cex=0.5, xlab="Effect size", ylab="Expected P value of Wilcoxon rank test")
plot(otu.p.all$diff.btw, otu.p.all$wi.ep, log="y", pch=19,
     main="Volcano", cex=0.5, xlab="Difference between the groups", ylab="Expected P value of Wilcoxon rank test")
abline(h=0.05, lty=2,lwd=3, col='red')
x_limits <- c(NA, -0.7)
ggplot(data=otu.p.all, aes(x=effect, y=-log10(wi.ep), color=p.effect)) + geom_point(size=1.5, alpha=0.5) +
  scale_x_continuous(name='Effect size', limits=c(-0.85, 0.85)) +
  scale_y_continuous(name='Uncorrected P values (-log10)', limits=c(0, 6), breaks=c(seq(0,6,1))) +
  scale_color_manual(name='Abundence changes', labels=c('Decreased', 'Unchanged', 'Increased'), values=c('green', 'black','red')) +
  #geom_text(aes(label=rownames(otu.p.all)), check_overlap=TRUE) + 
  geom_text_repel(aes(label=sig.otu), color='black', xlim=x_limits) +
  geom_vline(xintercept=c(-0.5, 0.5), linetype=2, size=1, color='grey') +
  theme(plot.margin=unit(c(3,15,5,5), "mm"), 
        panel.background=element_rect(fill="white"),
        axis.line=element_line(color='black'),
        legend.position=c(0.95, 0.85), 
        legend.title=element_text(color='black', size=13), 
        legend.text=element_text(color='black', size=13), 
        axis.ticks=element_line(color='black', size=1), 
        axis.text=element_text(color='black', size=11),
        axis.title.x=element_text(size=15, vjust=-1), 
        axis.title.y=element_text(size=15, vjust=3))


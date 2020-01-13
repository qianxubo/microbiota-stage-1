library(randomForest)
library(ggplot2)

###Load and pre-progress the data
mydata <- read.csv('Table_S9.csv', skip=1, header=T)
mydata[1:5, 1:10]
genera <- mydata[,-c(1,2,3,4,6)]  #Delete taxonomy except Genus
genera_sum <- aggregate(genera[ ,-1], by=list(genera[ ,1]), sum)
rownames(genera_sum) <- genera_sum[, 1]
otu_g <- genera_sum[,-1]
otu.t <- data.frame(t(otu_g))
g.perc <- otu.t/rowSums(otu.t)
df <- round(g.perc*1000000, digits=0)
rowSums(df)

#Grouping factor
group_letter <- substring(rownames(df), 1, 1)
g <- ifelse(group_letter=="P", "JIA", "Control")
df$aa <- g
df <- df[,order(colnames(df), decreasing=F)]
str(df)
df$aa <- as.factor(df$aa)
str(df)

#var_error analysis
set.seed(69)
ww <- cbind(df[,-1], matrix(runif(100*nrow(df)), nrow(df), 100))
rr <- rfcv(ww, df[,1], cv.fold=10)
aa <- data.frame(rr$n.var, rr$error.cv*100)
aa  #Check OOB as the genera numbers increase

#Visualization
with(rr, plot(n.var, error.cv, log="x", type="o", lwd=2))
ggplot(data=aa, aes(x=rr.n.var, y=rr.error.cv...100)) + geom_line(size=0.5, color="black") +
  geom_vline(xintercept=c(12), linetype=2, size=1, color='grey') +
  scale_y_continuous(name='Error rates (%)', limits=c(32, 42), breaks=c(seq(32,42,2))) +
  scale_x_continuous(name='Variable numbers', limits=c(0, 100), breaks=c(seq(0,100,20))) +
  guides(color=FALSE) +
  theme(axis.title.x=element_text(size=15),
        axis.text.x=element_text(color='black', size=11),
        axis.title.y=element_text(size=15),
        axis.text.y=element_text(color='black', size=11),
        axis.line=element_line(),
        axis.ticks=element_line(color='black', size=1),
        panel.background=element_rect(fill="white"))

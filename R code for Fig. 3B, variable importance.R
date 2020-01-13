library(randomForest)
library(ggplot2) 

###Variable importance of genera
mydata <- read.csv('Table_S9.csv', skip=1, header=T)
mydata[1:5, 1:10]
genera <- mydata[,-c(1,2,3,4,6)]  #Delete taxonomy except Genus
genera[1:5, 1:10]
genera_sum <- aggregate(genera[ ,-1], by=list(genera[ ,1]), sum)
rownames(genera_sum) <- genera_sum[, 1]
otu_g <- genera_sum[,-1]
otu.t <- data.frame(t(otu_g))
df <- otu.t/rowSums(otu.t)
rowSums(df)
group_letter <- substring(rownames(df), 1, 1)
g <- ifelse(group_letter=="P", "JIA", "Control")
df$aa <- g
df <- df[,order(colnames(df), decreasing=F)] 
str(df)
df$aa <- as.factor(df$aa) 
str(df)
df[1:5, 1:10]
rowSums(df[,-1])
colnames(df)[1] <- "Groups"

set.seed(100) 
m.try <- tuneRF(df[-1], df$Groups, ntreeTry=500, stepFactor=1.5, improve=0.01, trace=TRUE)
m.try[m.try[,2] == min(m.try[,2]), 1]
set.seed(95)
rf <- randomForest(Groups~., data=df, ntree=500, mtry=9, na.action=na.omit, importance=TRUE)
rf 
im <- as.data.frame(rf$importance)
im <- im[order(im$MeanDecreaseGini, decreasing=T),]
im <- im[1:12,]
vi <- data.frame(row.names(im), im[,4]) 
colnames(vi) <- c("Genera", "Gini")
vi$Genera <- reorder(vi$Genera, vi$Gini)
ggplot(data=vi, aes(x=Genera, y=Gini, color=Genera)) + geom_point(size=4) +
  guides(color=FALSE) +
  scale_x_discrete(limits=rev(levels(vi$Genera))) +
  scale_y_continuous(name='Gini index') +
  theme(panel.background=element_rect(fill="white"), 
        plot.margin=unit(c(0.5,0.5,0.5,0.5), "cm"), 
        axis.title.y=element_text(size=15, vjust=5), 
        axis.title.x=element_text(size=15, vjust=-2), 
        axis.line=element_line(), 
        axis.text.x=element_text(color='black', size=13, angle=30, hjust=1, vjust=1),
        axis.text.y=element_text(color='black', size=13), 
        axis.ticks=element_line(color='black', size=1), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color="grey60", linetype="dashed"))
        

###Variable importance of OTUs
mydata <- read.csv('Table_S4.csv', skip=1, header=T)
mydata[1:5, 1:10]
myotu <- mydata[,-c(1:5)]  #Delete taxonomy except OTU
myotu[1:5, 1:10]
rownames(myotu) <- myotu[, 1]
otu <- t(myotu[,-1])
otu[1:5, 1:10]
d <- data.frame(otu/rowSums(otu))
rowSums(d)
d[1:5, 1:10]

group_letter <- substring(rownames(d), 1, 1)
g <- ifelse(group_letter=="P", "JIA", "Control")
str(g)
d$aa <- g 
d <- d[,order(colnames(d), decreasing=F)]
str(d$aa)
d$aa <- as.factor(d$aa)
str(d)
d[1:5, 1:5]
rowSums(d[,-1])
colnames(d)[1] <- "Groups"
d[1:5, 1:5]

set.seed(1000)
m.try <- tuneRF(d[-1], d$Groups, ntreeTry=500, stepFactor=1.5, improve=0.01, trace=TRUE)
m.try[m.try[,2] == min(m.try[,2]), 1]
set.seed(95)
rf <- randomForest(Groups~., data=d, ntree=500, mtry=25, na.action=na.omit, importance=TRUE)
rf
im <- as.data.frame(rf$importance)
im <- im[order(im$MeanDecreaseGini, decreasing=T),]
im <- im[1:8,]
vi <- data.frame(row.names(im), im[,4])

colnames(vi) <- c("OTUs", "Gini")
vi$OTUs <- reorder(vi$OTUs, vi$Gini)
ggplot(data=vi, aes(x=Gini, y=OTUs, color=OTUs)) + geom_point(size=4) +
  scale_x_continuous(name='Gini indices') +
  guides(color=FALSE) +
  theme(panel.background=element_rect(fill="white"), 
        plot.margin=unit(c(0.5,0.5,0.5,0.5), "cm"), 
        axis.title.y=element_text(size=15, vjust=5), 
        axis.title.x=element_text(size=15, vjust=-2), 
        axis.line=element_line(),  
        axis.text.x=element_text(color='black', size=13),
        axis.text.y=element_text(color='black', size=13),
        axis.ticks=element_line(color='black', size=1), 
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color="grey60", linetype="dashed"))


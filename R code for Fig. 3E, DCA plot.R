#Please load “R code for Figure 3E, DCA function.r” and and run "dca" function first. The last few steps can also be performed in Stata software.
library(randomForest)

###Load and pre-process the data
mydata <- read.csv('Table_S9.csv', skip=1, header=T)
mydata[1:5, 1:10]
genera <- mydata[,-c(1,2,3,4,6)]  #Delete taxonomy except Genus
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

###Random forest analysis, and then pick out 12 genera
set.seed(100)
m.try <- tuneRF(df[-1], df$Groups, ntreeTry=500, stepFactor=1.5, improve=0.01, trace=TRUE)
m.try[m.try[,2] == min(m.try[,2]), 1]
set.seed(95)
rf <- randomForest(Groups~., data=df, ntree=500, mtry=9, na.action=na.omit, importance=TRUE)
rf
im <- as.data.frame(rf$importance)
im <- im[order(im$MeanDecreaseGini, decreasing=T),]
im <- im[1:12,]
d <- df[,c("Groups", rownames(im))]


###DCA analysis
j=814
set.seed(j)
my.try <- tuneRF(d[-1], d$Groups, ntreeTry=500, stepFactor=1.5, improve=0.01, trace=F, plot=F)
m1 <- my.try[my.try[,2] == min(my.try[,2]), 1]
m2 <- m1[1]
set.seed(j)
myrf <- randomForest(Groups~., data=d, ntree=500, mtry=m2, na.action=na.omit, importance=TRUE)
myrf
mydca <- data.frame(myrf$votes[,2])
colnames(mydca)[1] <- "p"
group_letter <- substring(rownames(mydca), 1, 1)
mydca$group <- as.integer(ifelse(group_letter=="P", "1", "0"))
dca(data=mydca, outcome="group", predictors="p", smooth="TRUE", xstop=1)




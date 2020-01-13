library(randomForest)
library(ggplot2) 

mydata <- read.csv('Table_S9.csv', skip=1, header=T)
mydata[1:5, 1:10]
genera <- mydata[,-c(1,2,3,4,6)]  #Delete taxonomy except Genus
genera[1:5, 1:10]
unique(genera$Genus)  #94 genera
genera_sum <- aggregate(genera[ ,-1], by=list(genera[ ,1]), sum)
rownames(genera_sum) <- genera_sum[, 1]
df <- data.frame(t(genera_sum[,-1]))

group_letter <- substring(rownames(df), 1, 1)
g <- ifelse(group_letter=="P", "JIA", "Control")
df$aa <- g
df <- df[,order(colnames(df), decreasing=F)]
str(df)
df$aa <- as.factor(df$aa)
str(df)
df[1:5, 1:10]
colnames(df)[1] <- "Groups"
df[1:5, 1:10]
group_sum <- aggregate(df[ ,-1], by=list(df[ ,1]), sum)
rownames(group_sum) <- group_sum[, 1]
myvenn <- data.frame(t(group_sum[,-1]))
#write.csv(myvenn, file="venn.csv")

unique_in_JIA <- subset(myvenn, myvenn$Control < 1)
dim(unique_in_JIA)[1]
unique_in_Control <- subset(myvenn, myvenn$JIA < 1)
dim(unique_in_Control)[1]
dim(myvenn)[1] - dim(unique_in_JIA)[1] - dim(unique_in_Control)[1]  #Shared genera in the two groups

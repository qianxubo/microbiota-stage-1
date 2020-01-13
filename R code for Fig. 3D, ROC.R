library(randomForest)
library(ggplot2)
library(ROCR)

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
group_letter <- substring(rownames(df), 1, 1)  #Grouping factor
g <- ifelse(group_letter=="P", "JIA", "Control")
df$aa <- g
df <- df[,order(colnames(df), decreasing=F)]
str(df)
df$aa <- as.factor(df$aa)
str(df)
df[1:5, 1:10]
rowSums(df[,-1])
colnames(df)[1] <- "Groups"

###Select optimal mtry
set.seed(95)
m.try <- tuneRF(df[-1], df$Groups, ntreeTry=500, stepFactor=1.5, improve=0.01, trace=TRUE)
m.try[m.try[,2] == min(m.try[,2]), 1]

###Random forest analysis，and then pick out 12 genera
set.seed(95)
rf <- randomForest(Groups~., data=df, ntree=500, mtry=9, na.action=na.omit, importance=TRUE)
rf
im <- as.data.frame(rf$importance)
im <- im[order(im$MeanDecreaseGini, decreasing=T),]
im <- im[1:12,]
d <- df[,c("Groups", rownames(im))]

###Construct a random forest model using the 12 genera, and then extract vote probability for DCA analysis
set.seed(100)
my.try <- tuneRF(d[-1], d$Groups, ntreeTry=500, stepFactor=1.5, improve=0.01, trace=TRUE)
my.try[my.try[,2] == min(my.try[,2]), 1]
set.seed(100)
myrf <- randomForest(Groups~., data=d, ntree=500, mtry=2, na.action=na.omit, importance=TRUE)
myrf

###Plot ROC
K <- 10
CV <- function(n, Z=K, seed=100){
  hb1 <- rep(1:Z, ceiling(n/Z))[1:n]
  set.seed(seed)
  hb2 <- sample(hb1, n)
  mm <- list()
  for(i in 1:Z) mm[[i]] <- (1:n)[hb2 == i]
  return(mm)}
j = 274  #Set seeds
rn <- CV(nrow(d), K, j)
myroc <- matrix(NA, 10, 20)
myauc <- vector(mode="numeric", length=K)
for(i in 1:K){
    hm <- rn[[i]]
    test_set <- d[hm, ]
    train_set <- d[-hm, ]
    rf_roc <- randomForest(Groups~., data=train_set, ntree=500, na.action=na.omit, importance=TRUE)
    prob <- predict(rf_roc, test_set, type="prob")[,2]
    mypred <- prediction(prob, test_set$Groups)
    auc <- performance(mypred, measure="auc", x.measure="cutoff")
    myauc[i] <- auc@y.values[[1]]
    perf <- performance(mypred, "tpr", "fpr")
    myroc[1:length(perf@x.values[[1]]), 1:2 + 2*(i-1)] <- cbind(perf@x.values[[1]], perf@y.values[[1]])
  }
myroc
myauc
mean(myauc)
auc_se <- sd(myauc)/sqrt(10)


roc <- as.data.frame(apply(myroc, 2, function(x){x*100}))
roc$mx <- apply(roc[,c(1,3,5,7,9,11,13,15,17,19)], 1, mean)
roc$my <- apply(roc[,c(2,4,6,8,10,12,14,16,18,20)], 1, mean)
roc$se <- apply(roc[,c(2,4,6,8,10,12,14,16,18,20)], 1, sd)/sqrt(10)

ggplot(data=roc, aes(x=mx, y=my)) + geom_ribbon(aes(ymin=my-1.96*se, ymax=my+1.96*se), alpha=0.2) + geom_line() +
    annotate("segment", x=0, xend=100, y=0, yend=100, color="black") +
    geom_label(aes(x=79, y=4, label='AUC=79.75% \n 95% CI: 76.29%－83.20%'), color='black', size=5) +
    scale_x_continuous(name='False positive rate (%)', breaks=c(seq(0,100,20))) +
    scale_y_continuous(name='True positive rate (%)', breaks=c(seq(0,100,20))) +
    theme(axis.title.x=element_text(size=15),
          axis.text.x=element_text(color='black', size=11),
          axis.title.y=element_text(size=15),
          axis.text.y=element_text(color='black', size=11),
          axis.line=element_line(),
          axis.ticks=element_line(color='black', size=1),
          panel.background=element_rect(fill="white"))


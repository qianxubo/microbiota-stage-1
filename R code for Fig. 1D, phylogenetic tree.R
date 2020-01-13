library(ggplot2)
library(ggtree)
library(ggrepel)

rt=read.csv("Table_S6.csv", header=T, skip=1, row.names=1, comment.char = "")
dim(rt)
cls=list()
for(i in 1:nrow(rt)){
 otu=rownames(rt[i,])
 phylum=strsplit(as.character(rt$taxonomy[i]),"\\; |p\\_\\_")[[1]][3]
 cls[[phylum]]=c(cls[[phylum]], otu)
}
phylumNames=names(cls)
phylumNum=length(phylumNames)

tree <- read.tree("Tree_file.tree")
tree <- groupOTU(tree, cls)

ggtree(tree, layout='daylight', branch.length="none", aes(color=group)) +
  scale_color_brewer(name='Phyla', palette='Set2', direction=-1, breaks=phylumNames, labels=phylumNames)+ 
  theme(legend.position="right",
        legend.title=element_text(color='black', size=14),
        legend.text=element_text(color='black', size=13),
        legend.key.height=unit(8, 'mm'),
        legend.key.width=unit(8, 'mm')) + 
  geom_tiplab(aes(label=paste(" ", gsub("\\d+\\.\\d+|New\\.|Reference|CleanUp\\.", "", label), sep="")), size=3)


#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=T)
if(length(args) < 5)
{
  stop('Usage: freqsFragments.R <tree> <clusterList> <hmmOut> <outNameFull> <outNameFragments>')
}

suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(phangorn))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(readr))

trName <- args[1]
clsName <- args[2]
hmmName <- args[3] 
#fullName <- paste('../',args[4],sep='')
#fragName <- paste('../',args[5],sep='')
fullName <- args[4]
fragName <- args[5]

tr <- read.tree(trName)
tr <- midpoint(tr)
tl <- tr$tip.label
tg <- tl %>% strsplit(.,"_") %>% lapply(.,tail,1) %>% unlist %>% as.integer
cl <- read.table(clsName,header=F,row.names=NULL,sep="\t")
cl.tb <- as.data.frame.table(table(cl[,2]))
names(cl.tb) <- c("Cluster","Count")
annot <- cbind(Repr=tl[match(cl.tb[,1],tg)],cl.tb)
hmm <- read_lines(hmmName)

lhmm <- length(hmm)
ids <- rep("",lhmm)
clus <- rep(-1,lhmm)
for(i in 1:lhmm){
  l <- hmm[i]
  if(substr(l,1,1)!="#"){
    f <- strsplit(l,"\\s+")[[1]]
    ids[i] <- f[3]
    clus[i] <- strsplit(f[1],".",fixed=TRUE) %>% unlist %>% tail(.,1) %>% as.integer
  }
}
hmmdf <- data.frame(ids,clus)
if(max(hmmdf$clus,na.rm=T) != -1)
{
  hmmdf <- hmmdf[hmmdf[,1]!="",]
  hmmdf <- hmmdf[!duplicated(hmmdf[,1]),]
  hmmdf.tbl <- as.data.frame.table(table(hmmdf[,2]))
  hmmdf.tbl[,1] <- as.integer(hmmdf.tbl[,1])
  annot$Fragments <- 0
  annot$Fragments[match(hmmdf.tbl[,1],annot[,2])] <- hmmdf.tbl[,2]
  annot$Total <- annot$Count+annot$Fragments

  df <- fortify(tr)
  df$tipstats <- NA
  d1 <- df
  d2 <- df
  d2$tipstats[d2$isTip] <- annot[match(d2$label[d2$isTip],annot[,1]),]$Fragments
  d1$panel <- 'Tree'
  d2$panel <- 'FragmentCount'
  d1$panel <- factor(d1$panel, levels=c("Tree", "FragmentCount"))
  d2$panel <- factor(d2$panel, levels=c("Tree", "FragmentCount"))
  p <- ggplot(mapping=aes(x=x, y=y)) + facet_grid(.~panel, scale="free_x") + theme_tree2()
  ffdt <- p+geom_tree(data=d1) + geom_point(data=d2, aes(x=tipstats)) 
  ggsave(fragName,ffdt,device="svg",width=8,height=5,units="in")  
}

df <- fortify(tr)
df$tipstats <- NA
d1 <- df
d2 <- df
d2$tipstats[d2$isTip] <- annot[match(d2$label[d2$isTip],annot[,1]),]$Count
d1$panel <- 'Tree'
d2$panel <- 'FullCount'
d1$panel <- factor(d1$panel, levels=c("Tree", "FullCount"))
d2$panel <- factor(d2$panel, levels=c("Tree", "FullCount"))
p <- ggplot(mapping=aes(x=x, y=y)) + facet_grid(.~panel, scale="free_x") + theme_tree2()
fdt <- p+geom_tree(data=d1) + geom_point(data=d2, aes(x=tipstats)) 
ggsave(fullName,fdt,device="svg",width=8,height=5,units="in")


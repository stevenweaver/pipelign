#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=T)
if(length(args) < 3)
{
  stop('Usage: freqsLong.R <tree> <clusterList> <outname>')
}

suppressPackageStartupMessages(library(ape))

suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(phangorn))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(readr))

trName <- args[1]
clsName <- args[2]
#oName <- paste('../',args[3],sep='')
oName <- args[3]


tr <- read.tree(trName)
tr <- midpoint(tr)
tl <- tr$tip.label
tg <- tl %>% strsplit(.,"_") %>% lapply(.,tail,1) %>% unlist %>% as.integer
cl <- read.table(clsName,header=F,row.names=NULL,sep="\t")
cl.tb <- as.data.frame.table(table(cl[,2]))
names(cl.tb) <- c("Cluster","Count")
annot <- cbind(Repr=tl[match(cl.tb[,1],tg)],cl.tb)

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
dt <- p+geom_tree(data=d1) + geom_point(data=d2, aes(x=tipstats)) 
ggsave(oName,dt,device="svg",width=8,height=5,units="in")

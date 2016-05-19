#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=T)
if(length(args) < 2)
{
  stop('Usage: lengthDistribution.R <seqFile> <outName>')
}

suppressPackageStartupMessages(library(seqinr))

seqName <- args[1]
outName <- args[2]

seqs <- read.fasta(seqName)

numSeqs <- length(seqs)

seqLens <- rep(NA,numSeqs)

for(i in 1:numSeqs)
{
  seqLens[i] <- length(seqs[[i]])
}

svg(outName)
hist(seqLens,main='Sequence length distribution',xlab='Lengths',ylab='number of sequences')
invisible(dev.off())

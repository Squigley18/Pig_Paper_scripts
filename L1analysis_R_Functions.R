#!/bin/Rscript

#librarys used in R analysis

library(tidyverse) #Tidyverse library is needed for the use of functions such as ggplot
library(GenomicFeatures) #GenomicFeatures library is required to load the TXdb object for the segment graphs 
library(GenomicRanges) #GenomicRanges library is required to create GenomicRanges objects 
library(ggbio) #Ggbio library is required for annotating the genes within the segment graph 
library(biomaRt) #Biomart library is required for accessing the gene data of the pig X chromosome 
library(ggplot2)
library(svglite)

GRfromDF <- function(hit){
  makeGRangesFromDataFrame(hit,
                           keep.extra.columns=TRUE,
                           ignore.strand=FALSE,
                           seqinfo=NULL,
                           seqnames.field=c("seqnames"),
                           start.field="start",
                           end.field=c("end"),
                           strand.field="strand",
                           starts.in.df.are.0based=FALSE )  
}

overlaps.to.L1 <- function(tab){
  subject <-  tab
  colnames(subject) <- c("start","end")
  seqnames <- "X"
  strand <- "*" 
  hits <- cbind(subject,seqnames,strand)
  hitsGF <- GRfromDF(hits)
  L1.acc <- make.windows(1,7878)
  overlap.hits <- countOverlaps(L1.acc,hitsGF)
  overlaps <- as.data.frame(overlap.hits)
  xaxis <- seq(1:7878)
  overlap.df <- as.data.frame(cbind(overlap.hits,xaxis))
  overlap.df <- overlap.df %>% dplyr::rename(yaxis = overlap.hits)
  transposase <- as.data.frame(cbind(start=2324,stop=3214,yaxis=-5))
  reverse.transcriptase <- as.data.frame(cbind(start=3740,stop=7558,yaxis=-5))
  repeat.1 <- as.data.frame(cbind(start=1230,stop=2080,yaxis=-5))
  repeat.2 <- as.data.frame(cbind(start=400,stop=680,yaxis=-5))
  UTR5 <- as.data.frame(cbind(start=15,stop=2330,yaxis=-11))
  UTR3 <- as.data.frame(cbind(start=7595,stop=7865,yaxis=-11))
  LINE <- as.data.frame(cbind(start=15,stop=7865,yaxis=-11))
  
  overlap.hit.plot <- ggplot(overlap.df,aes(x=xaxis,na.rm=TRUE))+
    geom_line(aes(y=yaxis),na.rm=TRUE)+
    geom_rect(data=LINE,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="darkblue", col=NA, inherit.aes=F)+
    geom_rect(data=transposase,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="red", col=NA, inherit.aes=F)+
    geom_rect(data=reverse.transcriptase,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="red", col=NA, inherit.aes=F)+
    geom_rect(data=UTR5,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="grey39", col=NA, inherit.aes=F)+
    geom_rect(data=UTR3,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="grey39", col=NA, inherit.aes=F)+
    geom_rect(data=repeat.1,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="mediumblue", col=NA, inherit.aes=F)+
    geom_rect(data=repeat.2,aes(xmin=start,xmax=stop,ymin=yaxis,ymax=yaxis+3),fill="mediumblue", col=NA, inherit.aes=F)+
    scale_x_continuous(name="L1 accession region (bp)",label=c(0,2500,5000,7500), breaks=c(0,2500,5000,7500))+
    theme_classic()+
    ylab("Number of BLAST sequence similarities")+
    annotate("text", x = c(3940,7730,1172,540,1655,5649,2769), y = c(-7,-7,-7,-1,-1,-1,-1),
             label = c("L1", "3'UTR", "5'UTR","repeat","repeat","reverse transcriptase", "transposase") , color="black", size=3)+
    theme(axis.title.x = element_text(color = "grey20", size = 9.5, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.text.x=element_text(size=9.5),axis.text.y=element_text(size=9.5),
          axis.title.y = element_text(color = "grey20", size = 9.5, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          title = element_text(color = "grey17", size= 10, angle = 0, hjust = .5, vjust = 0, face = "bold"))
}

make.windows = function(window.size, chr.size){
  starts = seq(1, chr.size-window.size+1, window.size)
  ends   = seq(window.size, chr.size, window.size)
  GRanges(seqnames="X",
          ranges=IRanges(
            start= starts,
            end  = ends),
          seqlengths=c(X=chr.size))
}

manipulate.data.for.overlap.ideogram.plot <- function(mydf){
  hit.gr <- GRfromDF(mydf)
  window.3000 <- make.windows(3000,126e6)
  hits <- countOverlaps(window.3000,hit.gr)
  xaxis <- seq(1:42000)
  counts <- as.data.frame(cbind(xaxis,hits))
  counts$ideogram <- counts$hits
  counts$ideogram[counts$hits >= 1] <- counts$xaxis[counts$hits >= 1]
  counts$ideogram[counts$ideogram==0] <- NA
  counts
}


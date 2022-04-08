#!/bin/Rscript

args = commandArgs(trailingOnly=TRUE) #allow arguments to be passed to R script

source("/home/sq16564/L1analysis_R_Functions.R")

all.extended.hit <-  read.table("all.extended.positions", header = FALSE, sep = "",
                       numerals = c("allow.loss", "warn.loss","no.loss"),
                       stringsAsFactors = TRUE)
extended.starts <-  read.table("start.extended.positions", header = FALSE, sep = "",
                                numerals = c("allow.loss", "warn.loss","no.loss"),
                                stringsAsFactors = TRUE)
extended.ends <-  read.table("end.extended.positions", header = FALSE, sep = "",
                               numerals = c("allow.loss", "warn.loss","no.loss"),
                               stringsAsFactors = TRUE)

all.extended.hit <- all.extended.hit %>% dplyr::rename(hit.start = V1,hit.stop = V2,s.start = V3,s.stop = V4,q.start = V5,q.stop = V6) 
extended.starts <- extended.starts %>% dplyr::rename(hit.start = V1,hit.stop = V2,s.start = V3,s.stop = V4,q.start = V5,q.stop = V6) 
extended.ends <- extended.ends %>% dplyr::rename(hit.start = V1,hit.stop = V2,s.start = V3,s.stop = V4,q.start = V5,q.stop = V6) 

all.extend.in.L1 <- all.extended.hit[,c(3,4)]
all.extend.in.L1[all.extend.in.L1$s.start >= all.extend.in.L1$s.stop,c("s.start","s.stop")] <- all.extend.in.L1[all.extend.in.L1$s.start >= all.extend.in.L1$s.stop,c("s.stop","s.start")]
start.extend.in.L1 <- extended.starts[,c(3,4)]
start.extend.in.L1[start.extend.in.L1$s.start >= start.extend.in.L1$s.stop,c("s.start","s.stop")] <- start.extend.in.L1[start.extend.in.L1$s.start >= start.extend.in.L1$s.stop,c("s.stop","s.start")]
end.extend.in.L1 <- extended.ends[,c(3,4)]
end.extend.in.L1[end.extend.in.L1$s.start >= end.extend.in.L1$s.stop,c("s.start","s.stop")] <- end.extend.in.L1[end.extend.in.L1$s.start >= end.extend.in.L1$s.stop,c("s.stop","s.start")]

all.plot <- overlaps.to.L1(all.extend.in.L1)
ggsave("all.extended.hits.over.L1.svg",all.plot,unit="mm",height=85,width=120)
start.plot <- overlaps.to.L1(start.extend.in.L1)
ggsave("start.extended.hits.over.L1.svg",start.plot,unit="mm",height=85,width=120)
end.plot <- overlaps.to.L1(end.extend.in.L1)
ggsave("end.extended.hits.over.L1.svg",end.plot,unit="mm",height=85,width=120)

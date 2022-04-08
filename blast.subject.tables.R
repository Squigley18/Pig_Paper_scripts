#!/bin/Rscript

#librarys used in R analysis

library(tidyverse) #Tidyverse library is needed for the use of functions such as ggplot
library(GenomicFeatures) #GenomicFeatures library is required to load the TXdb object for the segment graphs 
library(GenomicRanges) #GenomicRanges library is required to create GenomicRanges objects 
library(ggbio) #Ggbio library is required for annotating the genes within the segment graph 
library(biomaRt) #Biomart library is required for accessing the gene data of the pig X chromosome 


args = commandArgs(trailingOnly=TRUE) #take arguments passed to script 


if (length(args)==0) { # if no arguments passed to script then kill script 
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) { #otherwise with arguments run the following functions


subject <-  read.table(args[1], header = FALSE, sep = "",
                            numerals = c("allow.loss", "warn.loss","no.loss"),
                            stringsAsFactors = TRUE) # read the table of blast subject and query start and stop positions into R with LASTZ hit start and stop positions

subject <- subject %>% dplyr::rename(hit.start = V1,hit.stop = V2,s.start = V3,s.stop = V4,q.start = V5,q.stop = V6) 

subject$xstart <- subject$hit.start + subject$q.start - 1   # calculate the start position of the subject sequence in the X chromosome and offset by one where the subject start position is calculated from position 1 in the LASTZ hit  
subject$xstop <- subject$hit.start + subject$q.stop - 1 # calculate the stop position of the subject sequence in the X chromosome and offset by one where the subject stop position is calculated from position 1 in the LASTZ hit

chromosome <- "X" 
x.subject <- as.data.frame(cbind(chromosome,subject$xstart,subject$xstop))
 names(x.subject) <- NULL #create table of subject start and stop positions in the X chromosome without a header to input into perl script 

write.table(subject,paste0(args[1],".blast.hits"),sep="\t",row.names=FALSE,quote=FALSE)


 write.table(x.subject,paste0(args[1],".chromosome.hits"),sep="\t",row.names = FALSE,quote = FALSE)
 
 }

library(rtracklayer)
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

repeats_10.2and11.1 <-  read.csv("NCBI.repeat.conversion.csv")
repeats_10.2and11.1 <- repeats_10.2and11.1[grep("chrX", repeats_10.2and11.1$mapped_id), ]
corrected.hits <- repeats_10.2and11.1[,c("mapped_id", "mapped_start", "mapped_stop", "source_start", "source_stop")]
colnames(corrected.hits) <- c("seqnames", "start", "end", "source_start", "source_stop")
corrected.hits$length <- corrected.hits$end - corrected.hits$start
full.hits <- corrected.hits[corrected.hits$length >= 5000, ]
full.L1.gr <- GRfromDF(full.hits)
start <- start(full.L1.gr)
end <- end(full.L1.gr)
og.start <- mcols(full.L1.gr)$source_start
og.end <- mcols(full.L1.gr)$source_stop
seqnames <- "chrX"
full.L1 <- as.data.frame(cbind(seqnames, start, end, og.start, og.end))
full.L1 <- transform(full.L1, end = as.integer(end), 
                     start = as.integer(start))
full.L1.gr <- GRfromDF(full.L1)

Rm99.L1s <- read.table("Rm99.L1.hits", header = FALSE)
colnames(Rm99.L1s) <- c("seqnames","start","end")
Rm95.L1s <- read.table("Rm95.L1.hits", header = FALSE)
colnames(Rm95.L1s) <- c("seqnames","start","end")
L1 <- rbind(Rm99.L1s,Rm95.L1s)
L1 <- unique(L1)
my.L1.gf <- GRfromDF(L1)
my.L1.gf <- reduce(my.L1.gf)
start <- start(my.L1.gf)
end<- end(my.L1.gf)
seqnames <- "chrX"
my.L1s <- as.data.frame(cbind(seqnames, start, end))
my.L1s <- transform(my.L1s, end = as.integer(end), 
                    start = as.integer(start))

my.L1s.gr <- GRfromDF(my.L1s)
my.start.L1.gr <- resize(my.L1s.gr,width=width(my.L1s.gr)+7878,fix='end')
my.end.L1.gr <- resize(my.L1s.gr,width=width(my.L1s.gr)+7878,fix='start')
extended.start <- as.data.frame(start(my.start.L1.gr))
extended.end <- as.data.frame(end(my.end.L1.gr))


find.start.overlaps <- findOverlaps(query = my.start.L1.gr, subject = full.L1.gr, type = 'any')

overlap.start.df <-  data.frame(extended.start[queryHits(find.start.overlaps),],my.L1s[queryHits(find.start.overlaps),], full.L1[subjectHits(find.start.overlaps),])

colnames(overlap.start.df) <- c("extended start","seqnames","my start","my end","seq2","11.1 start", "11.1 end", "10.2 start", "10.2 end")

find.end.overlaps <- findOverlaps(query = my.end.L1.gr, subject = full.L1.gr, type = 'any')

overlap.end.df <-  data.frame(extended.end[queryHits(find.end.overlaps),],my.L1s[queryHits(find.end.overlaps),], full.L1[subjectHits(find.end.overlaps),])

colnames(overlap.end.df) <- c("extended end","seqnames","my start","my end","seq2","11.1 start", "11.1 end", "10.2 start", "10.2 end")


no.apparent.overlap <- subset(my.L1s, !(start %in% overlap.start.df$"my start") & !(start %in% overlap.end.df$"my start"))
no.apparent.overlap$extended.start <- no.apparent.overlap$start - 7878
no.apparent.overlap$extended.end <- no.apparent.overlap$end + 7878


write.table(overlap.start.df,"overlap.start.L1.txt",sep="\t",
            row.names = FALSE,quote = FALSE)
write.table(overlap.end.df,"overlap.end.L1.txt",sep="\t",
            row.names = FALSE,quote = FALSE)
write.table(no.apparent.overlap,"no.overlap.L1.txt",sep="\t",
            row.names = FALSE,quote = FALSE)

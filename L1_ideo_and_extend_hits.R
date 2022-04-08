#!/bin/Rscript

args = commandArgs(trailingOnly=TRUE) #allow arguments to be passed to R script

source("/home/sq16564/L1analysis_R_Functions.R")

Rm99.hits <- read.table(args[1], header = FALSE)
colnames(Rm99.hits) <- c("seqnames","start","end")
Rm95.hits <- read.table(args[2], header = FALSE)
colnames(Rm95.hits) <- c("seqnames","start","end")
ideo.99.L1 <- manipulate.data.for.overlap.ideogram.plot(Rm99.hits)
ideo.95.L1 <- manipulate.data.for.overlap.ideogram.plot(Rm95.hits)

chromosome.ideogram <- data.frame(x = 1:42000, y = c(rep.int(1,20500),
                                                       rep.int(NA,800),rep.int(1,20700)),
                                    z = c(rep.int(2,20500),
                                          rep.int(NA,800),rep.int(2,20700)))
  
Ideogram <- ggplot(chromosome.ideogram, aes(x,y),na.rm = TRUE)+ 
  geom_path(size = 4, lineend = "round",colour="gray87")+
  geom_path(aes(x,z),size = 4, lineend = "round",colour="gray87")+
  geom_segment(data = ideo.95.L1,aes(x=ideogram,xend=ideogram+25,
                                     y=1,yend=1),size=3.5,colour="#1026EB")+
  geom_segment(data = ideo.99.L1,aes(x=ideogram,xend=ideogram+25,
                                     y=2,yend=2),size=3.5,colour="#0A0A0A")+
  scale_x_continuous(name="Location on chromosome (Mbp)", label=c(0,30,60,90,120), breaks=c(0,10000,20000,30000,40000))+
  theme_classic()+
  theme(axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

ggsave("ideogram.my.L1.svg",Ideogram,unit="mm",height=85,width=120)
		
Rm99.L1s <- read.table(args[1], header = FALSE)
colnames(Rm99.L1s) <- c("seqnames","start","end")
Rm95.L1s <- read.table(args[2], header = FALSE)
colnames(Rm95.L1s) <- c("seqnames","start","end")
L1 <- rbind(Rm99.L1s,Rm95.L1s)
L1 <- unique(L1)
my.L1.gf <- GRfromDF(L1)
my.L1.gf <- reduce(my.L1.gf)
start <- start(my.L1.gf)
end<- end(my.L1.gf)
seqnames <- "chrX"
my.L1s <- as.data.frame(cbind(seqnames, start, end))
my.L1s$start <- as.numeric(as.character(my.L1s$start))
my.L1s$end <- as.numeric(as.character(my.L1s$end))

my.L1s.gr <- GRfromDF(my.L1s)
my.start.L1.gr <- resize(my.L1s.gr,width=width(my.L1s.gr)+7878,fix='end')
my.end.L1.gr <- resize(my.L1s.gr,width=width(my.L1s.gr)+7878,fix='start')
extended.start <- as.data.frame(start(my.start.L1.gr))
extended.end <- as.data.frame(end(my.end.L1.gr))

chromosome <- "X"
all.hit.extensions <- as.data.frame(cbind(chromosome,my.L1s,extended.start,extended.end))
start.extension <- all.hit.extensions[,c(1,5,4)]
colnames(start.extension) <- c("chromosome","start","stop")
end.extension <- all.hit.extensions[,c(1,3,6)]
colnames(end.extension) <- c("chromosome","start","stop")
all.extension <- rbind(start.extension,end.extension)
names(start.extension) <- NULL
names(end.extension) <- NULL
names(all.extension) <- NULL

write.table(start.extension,"extended_start_hits",sep="\t",row.names = FALSE,quote = FALSE)
write.table(end.extension,"extended_end_hits",sep="\t",row.names = FALSE,quote = FALSE)
write.table(all.extension,"extended_all_hits",sep="\t",row.names = FALSE,quote = FALSE)

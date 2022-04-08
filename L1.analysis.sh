#!/bin/bash
# File: PlotWrapper.sh
#
#$ -S /bin/bash
#$ -q all.q
#$ -m be
#$ -M sq16564@essex.ac.uk

grep "retrotransposon L1" refined.blast/Rm99**.final > Rm99.L1
grep "retrotransposon L1" refined.blast/Rm95**.final > Rm95.L1
grep "retrotransposon L1" refined.blast/ID99**.final > ID99.L1
grep "retrotransposon L1" refined.blast/ID95**.final > ID95.L1

cut -f1,4,5,6,7 Rm99.L1 > Rm99.L1.positions
cut -f1,4,5,6,7 Rm95.L1 > Rm95.L1.positions
cut -f1,4,5,6,7 ID99.L1 > ID99.L1.positions
cut -f1,4,5,6,7 ID95.L1 > ID95.L1.positions

rm Rm99.L1
rm Rm95.L1
rm ID99.L1
rm ID95.L1

sed 's/:/\t/g' Rm99.L1.positions > Rm99.L1.tab
sed 's/:/\t/g' Rm95.L1.positions > Rm95.L1.tab
sed 's/:/\t/g' ID99.L1.positions > ID99.L1.tab
sed 's/:/\t/g' ID95.L1.positions > ID95.L1.tab

rm Rm99.L1.positions
rm Rm95.L1.positions
rm ID99.L1.positions
rm ID95.L1.positions

cut -f6,7,9,10,11,12 Rm99.L1.tab > Rm99.L1.tab2
cut -f6,7,9,10,11,12 Rm95.L1.tab > Rm95.L1.tab2
cut -f6,7,9,10,11,12 ID99.L1.tab > ID99.L1.tab2
cut -f6,7,9,10,11,12 ID95.L1.tab > ID95.L1.tab2

rm Rm99.L1.tab
rm Rm95.L1.tab
rm ID99.L1.tab
rm ID95.L1.tab

/home/sq16564/./blast.subject.tables.R Rm99.L1.tab2
/home/sq16564/./blast.subject.tables.R Rm95.L1.tab2
/home/sq16564/./blast.subject.tables.R ID99.L1.tab2
/home/sq16564/./blast.subject.tables.R ID95.L1.tab2

rm Rm99.L1.tab2
rm Rm95.L1.tab2
rm ID99.L1.tab2
rm ID95.L1.tab2

mv Rm99.L1.tab2.chromosome.hits Rm99.L1.hits
mv Rm95.L1.tab2.chromosome.hits Rm95.L1.hits
mv ID99.L1.tab2.chromosome.hits ID99.L1.hits
mv ID95.L1.tab2.chromosome.hits ID95.L1.hits

perl /home/sq16564/perl.extract.tsv.files.pl Rm99.L1.hits
perl /home/sq16564/perl.extract.tsv.files.pl Rm95.L1.hits
perl /home/sq16564/perl.extract.tsv.files.pl ID99.L1.hits
perl /home/sq16564/perl.extract.tsv.files.pl ID95.L1.hits


/home/sq16564/L1_ideo_and_extend_hits.R Rm99.L1.hits Rm95.L1.hits

perl /home/sq16564/perl.extract.tsv.files.pl extended_all_hits
perl /home/sq16564/perl.extract.tsv.files.pl extended_start_hits
perl /home/sq16564/perl.extract.tsv.files.pl extended_end_hits

blastn  -query extended_all_hits.out.tsv -subject L1.sequence -out all.extended.blast -outfmt "6 qseqid sseqid sstart send qstart qend"
blastn  -query extended_start_hits.out.tsv -subject L1.sequence  -out start.extended.blast -outfmt "6 qseqid sseqid sstart send qstart qend"
blastn  -query extended_end_hits.out.tsv -subject L1.sequence  -out end.extended.blast -outfmt "6 qseqid sseqid sstart send qstart qend"


sed 's/:/\t/g' all.extended.blast > all.hit.positions
cut -f5,6,9,10,11,12 all.hit.positions > all.extended.positions
sed 's/:/\t/g' start.extended.blast > start.hit.positions
cut -f5,6,9,10,11,12 start.hit.positions > start.extended.positions
sed 's/:/\t/g' end.extended.blast > end.hit.positions
cut -f5,6,9,10,11,12 end.hit.positions > end.extended.positions


/home/sq16564/plot_extended_hits_over_L1.R


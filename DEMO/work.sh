#1.Mining single target gene from genome skimming data.
geneminer.py -1 skimming_data1.fq.gz -2 skimming_data2.fq.gz -rtfa shallow_ref/ITS.fasta -o ITS_out


2.Mining multiple target genes from genome skimming data.
#for FASTA format
geneminer.py -1 skimming_data1.fq.gz -2 skimming_data2.fq.gz -rtfa shallow_ref/ -o out1 
#for GenBank format
geneminer.py -1 skimming_data1.fq.gz -2 skimming_data2.fq.gz -rtgb shallow_ref.gb -o out2

#3.Mining multiple target genes from genome skimming data and evaluating  the assembly results.
geneminer.py -1 skimming_data1.fq.gz  -2 skimming_data2.fq.gz -rtfa shallow_ref/ -min 300 -max 5000 -limit_count auto -t 4 -bn 20 -o out3

#4.Mining Angiosperms353 genes from transcriptome data
geneminer.py -1 Arabidopsis_thaliana_sim1.fq.gz -2 Arabidopsis_thaliana_sim2.fq.gzz -rtfa ref_Arabidopsis_353 -k1 29 -k2 41 -t 4 -limit_count 2 -o Angiosperm353
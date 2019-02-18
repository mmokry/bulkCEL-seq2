#!/user/bin/perl -w
use strict;


#arguments 1 and 2 are R1 and R2 fastq files
open BC, $ARGV[2] or die;  #primers file



system ("gunzip -v *fastq.gz");


system ("cat *L001_R1_001.fastq *L002_R1_001.fastq *L003_R1_001.fastq *L004_R1_001.fastq > $ARGV[0].all");
system ("cat *L001_R2_001.fastq *L002_R2_001.fastq *L003_R2_001.fastq *L004_R2_001.fastq > $ARGV[1].all");



system ("perl 01_extract_barcodes_cell_seq.pl $ARGV[0].all $ARGV[1].all");


while (<BC>){
    my($name,$bc)=split(/\t/,$_);
    system ("./bwa aln -B 6 -q 0 -n 0.00 -k 2 -l 200 -t 6 HG19_genes.fa $name.R1.fastq > $name.R1.sai");
    system ("./bwa aln -B 0 -q 0 -n 0.04 -k 2 -l 200 -t 6 HG19_genes.fa $name.R2.fastq > $name.R2.sai");
    system ("./bwa sampe -n 100 -N 100 HG19_genes.fa $name.R1.sai $name.R2.sai $name.R1.fastq $name.R2.fastq > $name.sam");
    system (" perl 02_parse_sam_cell_seq.pl $name.sam ");

 



    }


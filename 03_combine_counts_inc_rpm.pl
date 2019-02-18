#!/user/bin/perl -w
use strict;


open OUT, '>AS_combined_raw_counts.txt' or die;
open OUTR, '>AS_combined_RPM_counts.txt' or die;


my%gene;
my%gener;
my%depth;
my $header = "gene\t";
foreach my $sample (@ARGV){
    open GENE, $sample;
    <GENE>;
    my $c=0;
    $depth{$sample}=0;
    while (my$line = <GENE>){
	chomp($line);
	my @r = split(/\t/,$line);
	$c++;
	$gene{$r[0]}=$r[0];
	$gener{$r[0]}=$r[0];
	if ($r[1]>4095){$r[1]=4095}
	$depth{$sample}+=-4096*(log(1-($r[1]/4096)));
    }
    print "$c genes in $sample.\n";
    print "$depth{$sample} mapped tags in $sample.\n\n";
    close(GENE);
}

print scalar(keys%gene)." genes in total\n";



foreach my $sample (@ARGV){
    $header .= $sample."\t";
    print "processing gene counts in $sample.\n";
    open GENE, $sample;
    <GENE>;
    my %gene_s;
    my %gener_s;
    
    while (my$line = <GENE>){
	chomp($line);
	my @r = split(/\t/,$line);
	if ($r[1]>4095){$r[1]=4095}
	$gene_s{$r[0]}=sprintf("%.0f"  , (-4096*(log(1-($r[1]/4096)))));
	$gener_s{$r[0]}=((-4096*(log(1-($r[1]/4096))))/$depth{$sample})*1000000;
    }
    foreach my $g (sort keys%gene){
        if ($gene_s{$g}){
    	    $gene{$g}.= "\t$gene_s{$g}";
    	    $gener{$g}.= "\t$gener_s{$g}";
    	    
	}
	else {
	    $gene{$g}.= "\t0";
	    $gener{$g}.= "\t0";
	}
    }
}





chomp($header);
print OUT $header ."\n";
print OUTR $header ."\n";

foreach my $g (sort keys%gene){
    print OUT "$gene{$g}\n";
    print OUTR "$gener{$g}\n";
    
}
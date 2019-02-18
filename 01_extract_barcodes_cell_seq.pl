#!/user/bin/perl -w
use strict;

open BC, 'primers.txt' or die;




while (my $line = <BC>){
open R1, $ARGV[0] or die;
open R2, $ARGV[1] or die;
    my$c_bc = 0;
    my %list;
    chomp($line);
    my @r=split(/\t/,$line);
    my $sample = $r[0];
    my $bc = $r[1];
    print "processing $line\n";
    open OUT, ">$sample".'.R2.fastq';
    open OUT2, ">$sample".'.R1.fastq';
    while (my$r1=<R1>){
        my $read = <R1>;
	my$r3=<R1>;
	my$r4=<R1>;
	my $umc = substr $read, 0,6;
	my $cel_bc = substr $read, 6,6;
	my$l1=<R2>;
	my$l2=<R2>;
	my$l3=<R2>;
	my$l4=<R2>;
	
	$l2 = substr $l2,0,32;
	$l4 = substr $l4,0,32;
	$l2.="\n";
	$l4.="\n";
	
	next unless ($cel_bc eq $bc);
	$c_bc++;
	my$id="$umc".'.'."$l2";
	print OUT "$l1$l2$l3$l4";
	print OUT2 "$r1$read$r3$r4";
	
    }
    print "all $c_bc\n";
}

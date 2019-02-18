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




























=head
open STAT, ">$ARGV[0]".'_mapping_length_'."$a_length".'.stat';
print STAT "barcode_primer\tname\treads_mapped\n";

my $time = time;
my $current = time;
my $elapsed = $current-$time;
$time = time;
print "$elapsed seconds since last timepoint\n";
############    extract the sequence after primer ################
while (my$seed = <BC>){
    print "processing $seed";
    chomp($seed);
    my @r = split(/\t/,$seed);
    $seed = $r[0];    
    open SEQ, $ARGV[0] or die;
    open OUT, ">$seed$ARGV[0]".'.seq' or die;
    my $count = 0;
	while (my $line = <SEQ>){
	    chomp($line);
	    my($n,$nn,$nnn,$nnnn,$nnnnn,$seq) = split(/\t/,$line);
	    next unless ($seq =~ s/^$seed//);
	    $seq = substr $seq, 0,$a_length;
    	    $count++;
    	    print OUT $seq."\n";	    
    }
    print "$count reads found for $r[1]\n";
    print STAT "$r[0]\t$r[1]\t$count\n";
    $current = time;
    $elapsed = $current-$time;
    $time = time;
    print "$elapsed seconds since last timepoint\n\n";
}
close STAT;

################    map the restriction sites and reads        ###########################


open BC, 'primers.txt' or die;
while (my$seed = <BC>){
    chomp($seed);
    my @r = split(/\t/,$seed);
    my $seed = $r[0];
    my $name = $r[1];
    print "extracting $seed\n";
######
    open CHRF, 'all_chromosomes.cod' or die 'cannot open chromosome length file cod format';
    my$chrN='';
    my %keychr;
    my %keyseq;
    my %keyseqR;
    my %keyseqRkey;
    my %keypos;
    while (<CHRF>){
	my@s=split/\t/;
	$chrN=$s[1];
	open CHR, "$chrN".'.fa' or die;
	my $seq;
	<CHR>;
	while (<CHR>){
	    chomp;
	    $seq.=$_;
	}
	print length($seq)."\n";
	$seq =~ tr/[a-z]/[A-Z]/; 
	while ($seq =~m/GATC/g){
	    my$key = substr $seq, pos $seq,$a_length;
	    $key =~ tr/[a-z]/[A-Z]/;
	    if ($keyseq{$key}){
		$keyseq{$key}=-100000000;
	    }
	    else {
		$keyseq{$key}=1;
	    }
	    if ($keyseqR{$key}){
		$keyseqR{$key}=-100000000;
	    }
	    my$keyR = substr $seq, ((pos $seq)-4-$a_length),$a_length;
	    $keyR = reverse($keyR);
	    $keyR =~ tr/ACGT/TGCA/;
	    $keyseqRkey{$key}=$keyR;
	    if ($keyseqR{$keyR}){
		$keyseqR{$keyR}=-100000000;
	    }
	    else {
		$keyseqR{$keyR}=1;
	    }
	    if ($keyseq{$keyR}){
		$keyseq{$keyR}=-100000000;
	    }
	    
	    $keypos{$key}=pos $seq;
	    $keychr{$key}=$chrN
	         
	}
    }
    $current = time;
    $elapsed = $current-$time;
    $time = time;
    print "$elapsed seconds since last timepoint\n";

#######################
#    open CHRF, 'all_chromosomes.cod' or die 'cannot open chromosome length file cod format';
#    $chrN='';
#    while (<CHRF>){
#	my@s=split/\t/;
#	$chrN=$s[1];
#	print "removing duplicates from $chrN\n";
#	open CHR, "$chrN".'.fa' or die;
#	my$seq2;
#	<CHR>;
#	while (<CHR>){
#	    chomp;
#	    $seq2.=$_;
#	}
#	my %testkeyseq;
#	$seq2 =~ tr/[a-z]/[A-Z]/; 
#	while ($seq2 =~m/GATC/g){
#	    my$key = substr $seq2, pos $seq2,$a_length;
#	    $key =~ tr/[a-z]/[A-Z]/;
#	    if ($keyseq{$key}){
#		$keyseq{$key}=-100000000;
#	    }
#	    else {
#		$keyseq{$key}=1;
#	    }
#	    if ($keyseqR{$key}){
#		$keyseqR{$key}=-100000000;
#	    }


#	    my$keyR = substr $seq2, ((pos $seq2)-4-$a_length),$a_length;
#	    $keyR = reverse($keyR);
#	    $keyR =~ tr/ACGT/TGCA/;
#	    $keyseqRkey{$key}=$keyR;
#	    if ($keyseqR{$keyR}){
#		$keyseqR{$keyR}=-100000000;
#	    }
#	    else {
#		$keyseqR{$keyR}=1;
#	    }
#	    if ($keyseq{$keyR}){
#		$keyseq{$keyR}=-100000000;
#	    }
#	}
#    }

    $current = time;
    $elapsed = $current-$time;
    $time = time;
    print "$elapsed seconds since last timepoint\n";




########################33
    open SEQ, "$seed$ARGV[0]".'.seq' or die;
    while (<SEQ>){
	chomp;
	if($keyseq{$_}){
	    $keyseq{$_}++;
	}
	if($keyseqR{$_}){
	    $keyseqR{$_}++;
	}
    }
    open OUT, ">$name".'_mapping_length_'."$a_length".'_'."$ARGV[0]".'.wig';
    print OUT "track type=wiggle_0 name=$name description=$name color=0,0,200 smoothingWindow=12 visibility=full autoScale=on windowingFunction=mean maxHeightPixels=40:40:20 color=0,0,200  priority=10\n";
    open OUTF, ">$name".'full_mapping_length_'."$a_length".'_'."$ARGV[0]".'.wig';
    print OUTF "track type=wiggle_0 name=$name description=$name color=0,0,200 smoothingWindow=12 visibility=full autoScale=on windowingFunction=mean maxHeightPixels=40:40:20 color=0,0,200  priority=10\n";
    
    
    open CHRF, 'all_chromosomes.cod' or die 'cannot open chromosome length file cod format';
    $chrN='';
    while (<CHRF>){
	my@s=split/\t/;
	$chrN=$s[1];


	print OUT "variableStep  chrom=$chrN\n";
	print OUTF "variableStep  chrom=$chrN\n";
	print "writing $chrN into wig file\n";

	foreach my $key (sort {$keypos{$a} <=> $keypos{$b}} keys %keypos){
	    next unless ($keychr{$key}eq$chrN);
#	    print "$keypos{$key}\t$keyseq{$key}\n";
	    if ($keyseqR{$keyseqRkey{$key}}>1){
		if ($keyseq{$key}<-10){
		    $keyseq{$key}=1
		}
		$keyseq{$key}--;
		$keyseq{$key}+=$keyseqR{$keyseqRkey{$key}};
	    }
	    if ($keyseq{$key} > 1){
		$keyseq{$key}-=1;
		print OUT "$keypos{$key}\t$keyseq{$key}\n";
		print OUTF "$keypos{$key}\t$keyseq{$key}\n";
#		print "$keypos{$key}\t$keyseq{$key}\n";
		delete $keypos{$key};
	    }
	    elsif ($keyseq{$key} == 1){
		print OUTF "$keypos{$key}\t0\n";
		delete $keypos{$key};
	    }
	    else {
		delete $keypos{$key}
	    }
	}
    $current = time;
    $elapsed = $current-$time;
    $time = time;
    print "$elapsed seconds since last timepoint\n";
    }
}


=cut
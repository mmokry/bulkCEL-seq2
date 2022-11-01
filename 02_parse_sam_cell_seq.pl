#!/user/bin/perl -w
use strict;


open IN, $ARGV[0] or die;
open OUT, ">$ARGV[0]".'.counts' or die;
open OUTUMI, ">$ARGV[0]".'.UMI' or die;
open OUTQC, ">$ARGV[0]".'.QC' or die;


my%gene;
my%check;
my%UMI;
my$reads=0;
my$mapped=0;
my$unique=0;


while (<IN>){
    chomp();
    my@r=split(/\t/,$_);
    next unless ($r[6]);
    $reads+=0.5;
    next unless ($_=~ m/ENST/);
    my@p=split(/ENST/,$r[2]);
    chop($p[0]);
    next if ($r[11]=~m/n/);
    next if ($r[11]=~m/c/);
    next if ($r[11]=~m/a/);
    next if ($r[11]=~m/t/);
    next if ($r[11]=~m/g/);
    next if ($r[11]=~m/N/);
#    next if ($r[11]=~m/n/);
#    next if ($r[11]=~m/n/);



    $UMI{$r[11]}++;
    $UMI{$r[11]}-=0.5;
    $mapped+=0.5;
    next if ($check{$p[0]}{$r[11]});
    $check{$p[0]}{$r[11]}=1;
    $gene{$p[0]}++;
    $unique++;
}


foreach my $g (sort keys%gene){
    print OUT "$g\t$gene{$g}\n"
}


foreach my $u (sort keys%UMI){
    print OUTUMI "$u\t$UMI{$u}\n"
}


my$mapped_perc=$mapped/$reads*100;
my$unique_perc=$unique/$mapped*100;


print OUTQC 'reads'."\t$reads\n";
print OUTQC 'mapped'."\t$mapped\n";
print OUTQC 'unique'."\t$unique\n";
print OUTQC 'mapped_perc'."\t$mapped_perc\n";
print OUTQC 'unique_perc'."\t$unique_perc\n";



#!/user/bin/perl -w
use strict;


open IN, $ARGV[0] or die;
open OUT, ">$ARGV[0]".'.counts' or die;
open OUTUMI, ">$ARGV[0]".'.UMI' or die;


my%gene;
my%check;
my%UMI;


while (<IN>){
    chomp();
    my@r=split(/\t/,$_);
    next unless ($r[6]);
    my@p=split(/ENST/,$r[2]);
    chop($p[0]);
    next if ($r[11]=~m/n/);
    $UMI{$r[11]}++;
    next if ($check{$p[0]}{$r[11]});
    $check{$p[0]}{$r[11]}=1;
    $gene{$p[0]}++;
}


foreach my $g (sort keys%gene){
    print OUT "$g\t$gene{$g}\n"
}


foreach my $u (sort keys%UMI){
    print OUTUMI "$u\t$UMI{$u}\n"
}
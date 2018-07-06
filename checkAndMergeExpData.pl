#!/usr/bin/env perl
#============================================================================
# Name        		: checkAndMergeExpData.pl
# Author      		: 
# Version     		: v1.00
# Created On  		: Fri May 20 16:11:53 2016
# Last Modified By	: 
# Last Modified On	: Fri May 20 16:11:53 2016
# Copyright   		: Copyright (C) 2016
# Description 		: 
#============================================================================

=pod

=head1 Usage

    perl checkAndMergeExpData.pl [option] <infile> <outfile>

    -h  display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,);
GetOptions("h"	=>\$opt_h);
if(@ARGV<2 || $opt_h) { usage(); }
my $infile=shift @ARGV;
my $outfile=shift @ARGV;

my @col_name=();
my @out_col=();
open $in,$infile or die "Cann't open file $infile ($!) \n";
while(<$in>)
{
    chomp;
    my $line=$_;
    if(/^\s*$/ || /^#/) { next; }
    my @F=split /\t/;
    my ($prj,$expFile,$sampleID)=@F[0,1,2];
    my $pExp=readOneFile($expFile,$sampleID);
    if(!@col_name){
        @col_name=@{$pExp->{'col_1'}};
        push @out_col,$pExp->{'col_1'};
        push @out_col,$pExp->{'col_2'};
    }else{
        my $isInconstant=0;
        for(my $i=0; $i < @{$pExp->{'col_1'}};$i++)
        {
            if($pExp->{'col_1'}->[$i] ne $col_name[$i])
            {
                printf "ERROR: $sampleID\t$expFile\n";
                $isInconstant=1;
                last
            }
        }
        if(!$isInconstant)
        {
            push @out_col,$pExp->{'col_2'};
        }
    }
}
my $ncol=@out_col;
my $nrow=@{$out_col[0]};
print "ncol:$ncol\n";
print "nrow:$nrow\n";

open $out,">",$outfile or die "!";
for(my $i=0;$i<$nrow;$i++)
{
    my @o=();
    for(my $j=0;$j<$ncol;$j++)
    {
        push @o,$out_col[$j]->[$i];
    }
    print $out join("\t",@o)."\n";
}

###################################################################

sub readOneFile
{
    my ($infile,$sampleID)=@_;
    my $in;
    open $in,$infile or die "Cann't open file $infile ($!) \n";
    my $i=0;
    my @col_1=();
    my @col_2=();
    while(<$in>)
    {
        chomp;
        my $line=$_;
        $i++;
        my @F=split /\t/;
        if($i==1){
            $F[2]=$sampleID;
        }else{
            $F[2]=$F[2] * 1e6;
        }
        push @col_1,$F[0];
        push @col_2,$F[2];
    }
    my $ret={"col_1"=>\@col_1,"col_2"=>\@col_2};
    return $ret;
}


sub usage
{
    die `pod2text $0`;
}

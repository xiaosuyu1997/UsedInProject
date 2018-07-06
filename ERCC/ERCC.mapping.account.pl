#!/usr/bin/env perl
#============================================================================
# Name        		: ERCC.mapping.account.pl
# Author      		: 
# Version     		: v1.00
# Created On  		: Mon Sep 28 22:30:40 2015
# Last Modified By	: 
# Last Modified On	: Mon Sep 28 22:30:40 2015
# Copyright   		: Copyright (C) 2015
# Description 		: 
#============================================================================

=pod

=head1 Usage

    perl ERCC.mapping.account.pl [option] 

    -s  sample [default: SAMPLE]
    -a  analyzed.bam [required]
    -b  nomapping.bam [required]
    -h  display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);

my ($in,$out);
my ($opt_h,$opt_s,$opt_a,$opt_b);
GetOptions("h"	=>\$opt_h,"s=s"=>\$opt_s,"a=s"=>\$opt_a,"b=s"=>\$opt_b);
if(@ARGV<0 || $opt_h) { usage(); }
if(!defined($opt_s)) { $opt_s="SAMPLE"; }
if(!defined($opt_a)) { usage(); }
if(!defined($opt_b)) { usage(); }
#my $outfile=shift @ARGV;
my $infile=shift @ARGV;

my $nSample=0;
my $nERCC=0;
my $nNomapping=0;
open $in,"samtools view $opt_a |" or die "$!";
while(<$in>)
{
    chomp;
    my $line=$_;
    if(/^\s*$/ || /^#/) { next; }
    my @F=split /\t/;
    my ($ID,$chr)=@F[0,2];
    if($chr=~/^ERCC/) { $nERCC++; }
    else{ $nSample++; }
}
open $in,"samtools view $opt_b |" or die "$!";
while(<$in>)
{
    chomp;
    my $line=$_;
    if(/^\s*$/ || /^#/) { next; }
    my @F=split /\t/;
    my ($ID,$chr)=@F[0,2];
    $nNomapping++;
}
my $nn=$nSample+$nNomapping;
print "sample\tnERCC\tnSample\tnNomapping\tnERCC:nSample\tpSample\tpNomapping\n";
print join("\t",$opt_s,$nERCC,$nSample,$nNomapping,sprintf("%4.6f",$nERCC/$nSample),sprintf("%4.2f",100*$nSample/$nn),sprintf("%4.2f",100*$nNomapping/$nn))."\n";

###################################################################




sub usage
{
    die `pod2text $0`;
}

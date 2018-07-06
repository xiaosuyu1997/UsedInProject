#!/usr/bin/env perl
#============================================================================
# Name        		: TCR.tracer.filtered_TCRs.table.pl
# Author      		: 
# Version     		: v1.00
# Created On  		: Wed Nov  4 18:35:02 2015
# Last Modified By	: 
# Last Modified On	: Wed Nov  4 18:35:02 2015
# Copyright   		: Copyright (C) 2015
# Description 		: 
#============================================================================

=pod

=head1 Usage

    perl TCR.tracer.filtered_TCRs.table.pl [option] <infile>

    -h  display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use Data::Dumper;

my ($in,$out);
my ($opt_h,);
GetOptions("h"	=>\$opt_h);
if(@ARGV<0 || $opt_h) { usage(); }
#my $outfile=shift @ARGV;
my $infile=shift @ARGV;


if(defined($infile))
{
	if($infile=~/\.gz$/) { open $in,"gzip -cd $infile |" or die "Cann't open file $infile ($!) \n"; }
	elsif($infile=~/\.bam$/) { open $in,"samtools view $infile |" or die "Cann't open file $infile ($!) \n"; }
	else{ open $in,$infile or die "Cann't open file $infile ($!) \n"; }
}else
{
	open $in,"-" or die "$!";
}
my $sampleID="SAMPLE";
my $TCRA_recomb="";
my $TCRB_recomb="";
## -1; "TCRA"; "TCRB"
my $detail_section=-1;
my %res=();
while(<$in>)
{
    chomp;
    my $line=$_;
    my @F=split /\t/;
    ### get the sample id
    if(/^----/)
    {
        $sampleID=<$in>;
        chomp $sampleID;
        <$in>;
        next;
    }
    ### total 
    if(/^TCRA recombinants: (.+)/)
    {
        $TCRA_recomb=$1;
        next;
    }
    if(/^TCRB recombinants: (.+)/)
    {
        $TCRB_recomb=$1;
        next;
    }
    if(/^#TCRA#/) { $detail_section="TCRA"; next; }
    if(/^#TCRB#/) { $detail_section="TCRB"; next; }
    ## one recombination
    if(/^##(.+?)##/)
    {
        my $oneRecomb={'iID'=>$1,'VSeg'=>'','JSeg'=>'','ID'=>'','TPM'=>'NA','Productive'=>'NA','StopCodon'=>'NA','InFrame'=>'NA'};
        while(<$in>)
        {
            chomp;
            if(/^$/) { last; }
            if(/^V segment:\t(.+)/) { $oneRecomb->{"VSeg"}=$1; }
            if(/^J segment:\t(.+)/) { $oneRecomb->{"JSeg"}=$1; }
            if(/^ID:\t(.+)/) { $oneRecomb->{"ID"}=$1; }
            if(/^TPM:\t(.+)/) { $oneRecomb->{"TPM"}=$1; }
            if(/^Productive:\t(.+)/) { $oneRecomb->{"Productive"}=$1; }
            if(/^Stop codon:\t(.+)/) { $oneRecomb->{"StopCodon"}=$1; }
            if(/^In frame:\t(.+)/) { $oneRecomb->{"InFrame"}=$1; }
        }
        if(!$res{$detail_section}) { $res{$detail_section}=[]; }
        push @{$res{$detail_section}},$oneRecomb;
    }
}

#print STDERR Dumper(%res);
print "TCRChain\tSample\tRecombinants\tiID\tVSeg\tJSeg\tID\tTPM\tProductive\tStopCodon\tInFrame\n";
foreach my $p (@{$res{"TCRA"}})
{
    print join("\t","TCRA",$sampleID,$TCRA_recomb,$p->{'iID'}, $p->{'VSeg'}, $p->{'JSeg'}, $p->{'ID'}, $p->{'TPM'}, $p->{'Productive'}, $p->{'StopCodon'}, $p->{'InFrame'})."\n";
}
foreach my $p (@{$res{"TCRB"}})
{
    print join("\t","TCRB",$sampleID,$TCRB_recomb,$p->{'iID'}, $p->{'VSeg'}, $p->{'JSeg'}, $p->{'ID'}, $p->{'TPM'}, $p->{'Productive'}, $p->{'StopCodon'}, $p->{'InFrame'})."\n";
}
###################################################################




sub usage
{
    die `pod2text $0`;
}

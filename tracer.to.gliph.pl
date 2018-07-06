#!/usr/bin/env perl
#============================================================================
# Name        		: tracer.to.gliph.pl
# Author      		: zhenglt@gmail.com
# Version     		: v1.00
# Created On  		: Thu Dec 10 09:10:04 2015
# Last Modified By	: 
# Last Modified On	: Thu Dec 10 09:10:04 2015
# Copyright   		: Copyright (C) 2015
# Description 		: 
#============================================================================

=pod

=head1 Usage

    perl tracer.to.gliph.pl [option] <infile> <outfile>

    -s  sampleID [default: SAMPLE]
    -h  display this help and exit

=cut

###-a  TPM threshold of alpha [default: 10]
###-b  TPM threshold of beta [default: 15]

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use Data::Table;
use List::Util qw(max);

my ($in,$out);
my ($opt_h,$opt_s,$opt_r,$opt_a,$opt_b);
GetOptions("h"	=>\$opt_h,"s=s"=>\$opt_s,"r"=>\$opt_r,"a=i"=>\$opt_a,"b=i"=>\$opt_b);
if(@ARGV<2 || $opt_h) { usage(); }
if(!defined($opt_s)) { $opt_s="SAMPLE"; }
##if(!defined($opt_a)) { $opt_a=10; }
##if(!defined($opt_b)) { $opt_b=15; }

my $infile=shift @ARGV;
my $outfile=shift @ARGV;

run($infile,$outfile);

###################################################################



sub run
{
    my ($in,$out);
    my ($infile,$outfile)=@_;
    my %clonotype=();
    my %clonotypeCount=();
    my @inputData=();
    my %cellClonotype=();
	my %betaClone=();

    open $out,">",$outfile or die "Cann't open file $outfile ($!) \n";
	my $in_table = Data::Table::fromTSV($infile, 1);
	my $next = $in_table->iterator();
    my $cid=0;
	print $out "CDR3b\tTRBV\tTRBJ\tCDR3a\tTRAV\tTRAJ\tPatient\tCounts\n";
	#print $out join("\t",$in_table->header,"C_cross","C_cross_b","C_cross_npatients","Beta_npatients")."\n";
	#print join("\t","CDR3_Beta","ncells","npatients","patients")."\n";
	while (my $row = $next->())
	{
        my $cellID=$row->{"Cell_Name"};
		my $patient="PATIENT";
		if($row->{"patient"}){ $patient=$row->{"patient"}; }
        my @idAlphaBeta4=($row->{"Identifier(Alpha1)"},$row->{"Identifier(Alpha2)"},
			$row->{"Identifier(Beta1)"},$row->{"Identifier(Beta2)"});
        my @cdr3AlphaBeta4=($row->{"CDR3(Alpha1)"},$row->{"CDR3(Alpha2)"},$row->{"CDR3(Beta1)"},$row->{"CDR3(Beta2)"});
        my @productiveAlphaBeta4=($row->{"Productive(Alpha1)"},$row->{"Productive(Alpha2)"},
			$row->{"Productive(Beta1)"},$row->{"Productive(Beta2)"});
        my @TPMAlphaBeta4=($row->{"TPM(Alpha1)"},$row->{"TPM(Alpha2)"},$row->{"TPM(Beta1)"},$row->{"TPM(Beta2)"});
		##my @cdr3ntAlphaBeta4=($row->{""},$row->{},$row->{},$row->{});
		##my @VDJAlphaBeta4=($row->{},$row->{},$row->{},$row->{});
		##my @inFrameAlphaBeta4=($row->{},$row->{},$row->{},$row->{});
		##my @stopCodonAlphaBeta4=($row->{},$row->{},$row->{},$row->{});
		for(my $j=0;$j<2;$j++){
			if($idAlphaBeta4[$j] ne "NA" && $productiveAlphaBeta4[$j]=~/True/i && $cdr3AlphaBeta4[$j]!~/Couldn/){
				my ($vseg_a,$jseg_a)=$idAlphaBeta4[$j]=~/^(TRAV.+)_.+(TRAJ.+?)$/;
				$cdr3AlphaBeta4[$j]=~s/FG.G$/F/;
				for(my $i=0;$i<2;$i++){
					if($idAlphaBeta4[2+$i] ne "NA" && $productiveAlphaBeta4[2+$i]=~/True/i && $cdr3AlphaBeta4[2+$i]!~/Couldn/){
						my ($vseg_b,$jseg_b)=$idAlphaBeta4[2+$i]=~/^(TRBV.+)_.+(TRBJ.+?)$/;
						$cdr3AlphaBeta4[2+$i]=~s/FG.G$/F/;
						print $out "$cdr3AlphaBeta4[2+$i]\t$vseg_b\t$jseg_b\t$cdr3AlphaBeta4[$j]\t$vseg_a\t$jseg_a\t$patient\n";
					}
				}
			}
		}
    }
}

sub usage
{
    die `pod2text $0`;
}

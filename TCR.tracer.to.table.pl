#!/usr/bin/env perl
#============================================================================
# Name        		: TCR.tracer.consensus.clone.pl
# Author      		: 
# Version     		: v1.00
# Created On  		: Wed Aug 30 08:30:11 2017
# Last Modified By	: 
# Last Modified On	: Wed Aug 30 08:30:11 2017
# Copyright   		: Copyright (C) 2017
# Description 		: 
#============================================================================

=pod

=head1 Usage

    perl TCR.tracer.consensus.clone.pl [option] <infile>

    -i  cloneID column name [default: "cloneID"]
    -d  sequence in fasta file. id line with format: "^>(cellID)\|[AB]\|.+?\|(seqID)\sproductive=(.+?)\s.+cdr3=(.+?)\s"
    -h  display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use Data::Table;
use Data::Dumper;
use Bio::Seq;
use Bio::SeqIO;

my ($in,$out);
my ($opt_h,$opt_d,$opt_i);
GetOptions("h"	=>\$opt_h,"d=s"=>\$opt_d,"i=s"=>\$opt_i);
if(@ARGV<1 || $opt_h) { usage(); }
###my $outfile=shift @ARGV;
my $infile=shift @ARGV;
if(!defined($opt_i)){ $opt_i="C_strict"; }
my $idx_cloneID=$opt_i;

my $in_table = Data::Table::fromTSV($infile, 1);
my $next = $in_table->iterator();
my @h=$in_table->header;

my %D=();

if($opt_d=~/\.gz$/){
    open $in,"gzip -cd $opt_d |" or die "Cann't open file $opt_d ($!) \n";
}else{
    open $in,"$opt_d" or die "Cann't open file $opt_d ($!) \n";
}
my $inseq = Bio::SeqIO->new( -fh     => $in,
                             -format => "fasta");
while (my $seq = $inseq->next_seq) {
   if($seq->display_id()=~/^(.+?)\|([AB])\|.+?\|(.+)$/)
   {
       my ($cellID,$chain,$combID)=($1,$2,$3);
       if($seq->desc()=~/productive=(.+?)\s.+cdr3=(.+?)\s/)
       {
           my ($productive,$cdr3)=($1,$2);
           #print join("\t",$cellID,$chain,$combID,$productive,$cdr3)."\n";
           $D{$cellID}->{$combID}={'productive'=>$productive,'cdr3'=>$cdr3,'seq'=>$seq->seq()};
       }
   }
}                         

print join("\t",@h,("A1_productive","A2_productive","B1_productive","B2_productive"),("A1_cdr3","A2_cdr3","B1_cdr3","B2_cdr3"),("A1_seq","A2_seq","B1_seq","B2_seq"))."\n";
for(my $i = 0; $i < $in_table->nofRow; $i++)
{
    my @F=$in_table->row($i);
    my $row = $next->();
  
    my @ID_AB=($row->{"A1"},$row->{"A2"},$row->{"B1"},$row->{"B2"});
    my $cloneID=$row->{$idx_cloneID};
    my $cellID=$row->{"exampleCell"};
    my @productive_AB=("NA","NA","NA","NA");
    my @cdr3_AB=("NA","NA","NA","NA");
    my @seq_AB=("NA","NA","NA","NA");
    for(my $ii=0;$ii<@ID_AB;$ii++){
        if(exists($D{$cellID}->{$ID_AB[$ii]})){
            $productive_AB[$ii]=$D{$cellID}->{$ID_AB[$ii]}->{'productive'};
            $cdr3_AB[$ii]=$D{$cellID}->{$ID_AB[$ii]}->{'cdr3'};
            $seq_AB[$ii]=$D{$cellID}->{$ID_AB[$ii]}->{'seq'};
        }
    }
    print join("\t",@F,@productive_AB,@cdr3_AB,@seq_AB)."\n";
}



###################################################################


sub usage
{
    die `pod2text $0`;
}

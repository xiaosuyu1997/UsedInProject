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

    -i  cloneID column name [default: "C_strict"]
    -f  minimum fraction of cells of the same clone [default: 0.5]
    -k  NOT only productive [default: OFF]
    -h  display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use Data::Table;
use Data::Dumper;

my ($in,$out);
my ($opt_h,$opt_f,$opt_i,$opt_k);
GetOptions("h"	=>\$opt_h,"i=s"=>\$opt_i,"f=f"=>\$opt_f,"k"=>\$opt_k);
if(@ARGV<1 || $opt_h) { usage(); }
###my $outfile=shift @ARGV;
my $infile=shift @ARGV;
if(!defined($opt_f)){ $opt_f=0.5; }
if(!defined($opt_i)){ $opt_i="C_strict"; }
my $idx_cloneID=$opt_i;

my $in_table = Data::Table::fromTSV($infile, 1);
my $next = $in_table->iterator();
my @h=$in_table->header;

my %D=();
for(my $i = 0; $i < $in_table->nofRow; $i++)
{
    my @F=$in_table->row($i);
    my $row = $next->();
  
    my @ID_AB=($row->{"Identifier(Alpha1)"},$row->{"Identifier(Alpha2)"},$row->{"Identifier(Beta1)"},$row->{"Identifier(Beta2)"});
    my @productive_AB=($row->{"Productive(Alpha1)"},$row->{"Productive(Alpha2)"},$row->{"Productive(Beta1)"},$row->{"Productive(Beta2)"});
    my $cloneID=$row->{$idx_cloneID};
    my $cellID=$row->{"Cell_Name"};
    if(!exists($D{$cloneID})){
        $D{$cloneID}={'N'=>0,'TRA'=>{},'TRB'=>{}};
    }
    $D{$cloneID}->{'N'}++;
    for(my $ii=0;$ii<2;$ii++)
    {
        if($ID_AB[$ii] ne "NA" && ( (!$opt_k && $productive_AB[$ii] eq "True") || $opt_k ) ){
            push @{$D{$cloneID}->{'TRA'}->{$ID_AB[$ii]}},$cellID;
        }
    }
    for(my $ii=2;$ii<4;$ii++)
    {
        if($ID_AB[$ii] ne "NA" && ( (!$opt_k && $productive_AB[$ii] eq "True") || $opt_k ) ){
            push @{$D{$cloneID}->{'TRB'}->{$ID_AB[$ii]}},$cellID;
        }
    }
    ### for searching a typical cell
    for(my $ii=0;$ii<4;$ii++)
    {
        #        if($ID_AB[$ii] ne "NA"){
        #    $D{$cloneID}->{'cell'}->{$cellID}->{$ID_AB[$ii]}=1;
        #    }
        $D{$cloneID}->{'cell'}->{$cellID}->{$ID_AB[$ii]}++;
    }
}

print "cloneID\texampleCell\tN\tA1\tA2\tB1\tB2\tA1_n\tA2_n\tB1_n\tB2_n\n";
foreach my $cloneID (sort keys %D){
    my $N=$D{$cloneID}->{'N'};
    ### A1, A2, B1, B2
    my %consTRA=();
    my %consTRB=();
    my $ii=0;
    foreach (keys %{$D{$cloneID}->{'TRA'}}){
        my $seqCount=@{$D{$cloneID}->{'TRA'}->{$_}};
        if($seqCount > $opt_f*$N){
            $consTRA{$_}=$seqCount;
        }
        $ii++;
    }
    $ii=2;
    foreach (keys %{$D{$cloneID}->{'TRB'}}){
        my $seqCount=@{$D{$cloneID}->{'TRB'}->{$_}};
        if($seqCount > $opt_f*$N){
            $consTRB{$_}=$seqCount;
        }
        $ii++;
    }
    my @consSeq=();
    my @consSeqCount=();
    my @kk=keys %consTRA;
    if(@kk==0){
        push @consSeq,"NA","NA";
        push @consSeqCount,0,0;
    }elsif(@kk==1){
        push @consSeq,$kk[0],"NA";
        push @consSeqCount,$consTRA{$kk[0]},0;
    }elsif(@kk==2){
        push @consSeq,$kk[0],$kk[1];
        push @consSeqCount,$consTRA{$kk[0]},$consTRA{$kk[1]};
    }
    @kk=keys %consTRB;
    if(@kk==0){
        push @consSeq,"NA","NA";
        push @consSeqCount,0,0;
    }elsif(@kk==1){
        push @consSeq,$kk[0],"NA";
        push @consSeqCount,$consTRB{$kk[0]},0;
    }elsif(@kk==2){
        push @consSeq,$kk[0],$kk[1];
        push @consSeqCount,$consTRB{$kk[0]},$consTRB{$kk[1]};
    }
    ### example cell
    my %consSeq_l=();
    foreach (@consSeq) { $consSeq_l{$_}++; }
    my $eCell="NA";
    foreach my $cellID (keys %{$D{$cloneID}->{'cell'}})
    {
        #if($cellID eq "NTC209-20161228"){
        #    print STDERR Dumper($D{$cloneID}->{'cell'}->{$cellID}),"\n";
        #}
        my $s=0;
        my $t=0;
        foreach (keys %consSeq_l)
        {
            $t++;
            if(exists($D{$cloneID}->{'cell'}->{$cellID}->{$_}) && $D{$cloneID}->{'cell'}->{$cellID}->{$_}==$consSeq_l{$_}){
                $s++;
            }
        }
        if($s==$t && $s>0){
            $eCell=$cellID;
            last;
        }
    }
    print join("\t",$cloneID,$eCell,$N,@consSeq,@consSeqCount)."\n";
}

###################################################################


sub usage
{
    die `pod2text $0`;
}

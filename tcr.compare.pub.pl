#!/usr/bin/env perl
#============================================================================
# Name        		: tcr.compare.pub.pl
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

    perl tcr.compare.pub.pl [option] <infile>

    -h  display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use Data::Table;
use Data::Dumper;

my ($in,$out);
my ($opt_h,);
GetOptions("h"	=>\$opt_h);
if(@ARGV<0 || $opt_h) { usage(); }
###my $outfile=shift @ARGV;
my $infile=shift @ARGV;

my %db_file=("Warren2011"=>"/DBS/DB_temp/zhangLab/repseq/TCR/Warren.TCRb2010/TCRBclonotypes/all.aaCDR3.uniq",
             "Li2016"=>"/DBS/DB_temp/zhangLab/repseq/TCR/Shirley.Liu/Supplementary.Table.1.cononical.seqOnly.txt");
my %db=();
foreach (keys %db_file){
    $db{$_}={};
    readList($db{$_},$db_file{$_});
}    
###print Dumper(\%db);
###exit;

my @dnames=keys %db;
my $in_table = Data::Table::fromTSV($infile, 1);
my $next = $in_table->iterator();
my @h=$in_table->header;
##print join("\t",@h,"CDR3(Beta1)(".join("/",@dnames).")","CDR3(Beta2)(".join("/",@dnames).")")."\n";
print join("\t",@h,"public_dataset(".join("/",@dnames).")")."\n";
my %stat=();
for(my $i = 0; $i < $in_table->nofRow; $i++)
{
    my @F=$in_table->row($i);
    my $row = $next->();
    
    ####my @cdr3Beta2=($row->{"CDR3(Beta1)"},$row->{"CDR3(Beta2)"});
    my @cdr3Beta2=($row->{"CDR3_Beta"});
    ###my @productiveBeta2=($row->{"Productive(Beta1)"},$row->{"Productive(Beta2)"});
    my @addField=();
    for(my $k=0;$k<@cdr3Beta2;$k++){
        my $aaCDR3=$cdr3Beta2[$k];
        ####if($aaCDR3!~/Couldn/ && $aaCDR3 ne "NA" && $productiveBeta2[$k]=~/True/i){
        if($aaCDR3!~/Couldn/ && $aaCDR3 ne "NA"){
            $aaCDR3=~s/FG.G$/F/;
            my @aaCDR3InDB_bool=();
            my $inPub=0;
            my $_part=0;
            for(my $j=0;$j<@dnames;$j++){
                if(exists($db{$dnames[$j]}->{$aaCDR3})){
                    push @aaCDR3InDB_bool,$db{$dnames[$j]}->{$aaCDR3};
                    #$stat{$dnames[$j]}++;
                    #$inPub=1;
                    $_part+=2**$j;
                }else{
                    push @aaCDR3InDB_bool,"No";
                }
            }
            $stat{"part$_part"}++;
            ##if($inPub==0){ $stat{"this_study"}++; }
            push @addField,join("/",@aaCDR3InDB_bool);
        }else{
            push @addField,join("/",("NA")x@dnames);
            $stat{"NA"}++;
        }
    }
    print join("\t",@F,@addField)."\n";
}
print STDERR join("\t",@dnames)."\n";
foreach (sort keys %stat){
    print STDERR "$_\t$stat{$_}\n";
}

###################################################################



sub readList
{
    my $in;
    my ($pList,$infile)=@_;
    open $in,$infile or die "Cann't open file $infile ($!) \n";
    while(<$in>)
    {
        chomp;
        my $line=$_;
        if(/^\s*$/ || /^#/) { next; }
        my @F=split /\t/;
        $pList->{$F[0]}=$F[1];
    }
}

sub usage
{
    die `pod2text $0`;
}

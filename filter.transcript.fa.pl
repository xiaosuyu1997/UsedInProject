#!/usr/bin/env perl
#============================================================================
# Name        		: filter.transcript.fa.pl
# Author      		: 
# Version     		: v1.00
# Created On  		: Sat May 20 09:33:26 2017
# Last Modified By	: 
# Last Modified On	: Sat May 20 09:33:26 2017
# Copyright   		: Copyright (C) 2017
# Description 		: 
#============================================================================

=pod

=head1 Usage

    perl filter.transcript.fa.pl [option] <infile>

    -h  display this help and exit

=cut

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use Bio::SeqIO;

my ($in,$out);
my ($opt_h,);
GetOptions("h"	=>\$opt_h);
if(@ARGV<0 || $opt_h) { usage(); }
#my $outfile=shift @ARGV;
my $infile=shift @ARGV;

###################################################################

my $seq_in = Bio::SeqIO->new(-fh => \*STDIN, -format => "fasta");
my $seq_out = Bio::SeqIO->new(-fh => \*STDOUT, -format => "fasta");
while (my $inseq = $seq_in->next_seq)
{
    my $str_desc=$inseq->desc();
    my ($chr,$beg,$end,$strand)=$str_desc=~/:GRCh38:(.+?):(.+?):(.+?):(.+?)(\s|$)/;
    my ($gid)=$str_desc=~/gene:(.+?)(\s|$)/;
    my $tid=$inseq->display_id;
    my ($gname)="NA";
    if($str_desc=~/gene_symbol:(.+?)(\s|$)/){ $gname=$1; }
    my $gene_biotype="NA";
    if($str_desc=~/gene_biotype:(.+?)(\s|$)/){ $gene_biotype=$1; }
    my $transcript_biotype="NA";
    if($str_desc=~/transcript_biotype:(.+?)(\s|$)/){ $transcript_biotype=$1; }
    ##if($gname=~/^(PDCD1)$/)
    if($chr!~/^CHR/)
    {
        print STDERR join("\t",$tid,$gid,$gname,$chr,$beg,$end,$strand,$gene_biotype,$transcript_biotype)."\n";
        $seq_out->write_seq($inseq);
    }
}

sub usage
{
    die `pod2text $0`;
}

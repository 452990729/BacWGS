#!/usr/bin/perl

=head1 Name

tRNAscan_to_gff3.pl  --  convert result of tRNAscan to gff3 format

=head1 Description

Convert the raw result of tRNAscan to gff3 format, the input files
include the position files, and(or) the structure files.

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Mender:huangqf,huangqf@genomic.org.cn
  Version: 1.1,  Date: 2008-4-28
  Note:

=head1 Usage

  $perl tRNAscan_to_gff3.pl <tRNA_file> [tRNA_structure_file]
  --prefix <str>   set the prefix tag of gene name
  --verbose        output running progress information to screen  
  --help           output help information to screen  

=head1 Exmple

perl ./tRNAscan_to_gff3.pl kaikoall-test-tRNA.fa.tRNA kaikoall-test-tRNA.fa.tRNA.structure --prefix Bm > kaikoall-test-tRNA.fa.tRNA.gf

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;


my ($Prefix,$Verbose,$Help);
GetOptions(
	"prefix:s"=>\$Prefix,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$Prefix = $Prefix."_" if(defined $Prefix);
die `pod2text $0` if (@ARGV == 0 || $Help);

my $tRNA_file = shift;
my $tRNA_structure_file = shift;

my %trna;


parse_tRNA($tRNA_file,\%trna) if(-f $tRNA_file);

parse_structure($tRNA_structure_file,\%trna) if(-f $tRNA_structure_file);

##print Dumper \%trna;

##output in gff3 format
my $output;
my $mark = "0001";
foreach my $chr (sort keys %trna) {
	my $chr_p = $trna{$chr};
	foreach my $trna (sort keys %$chr_p) {
		my $tRNA_p = $chr_p->{$trna};
		my $tRNA_id = $Prefix."tRNA_".$mark;
		$output .= "$chr\ttRNAscan-SE\ttRNA\t$tRNA_p->{start}\t$tRNA_p->{end}\t$tRNA_p->{score}\t$tRNA_p->{strand}\t\.\tID=$tRNA_id;Type=$tRNA_p->{type};Anti-codon=$tRNA_p->{anticodon};";
		$output .= "structure=\"$tRNA_p->{structure}\"; " if(defined $tRNA_p->{structure});
		$output .= "\n";
		$mark++;
	}
}

print "##gff-version 3\n$output";

####################################################
################### Sub Routines ###################
####################################################


sub parse_structure {
	my $file = shift;
	my $data = shift;

	$/="\n\n";
	open IN,$file || die "fail $file";
	while (<IN>) {
		my ($seq,$tRNA) = ($1,$1.'.'.$2) if(/^(\w+)\.(\w+)/);
		next if(!exists $trna{$seq} || !exists $trna{$seq}{$tRNA});
		my $structure = $1 if(/\nStr: (\S+)/);
		my $sequence = $1 if(/\nSeq: (\S+)/); ## the tRNA DNA sequence
		$trna{$seq}{$tRNA}{structure} = $structure;
		$trna{$seq}{$tRNA}{sequence} = $sequence;
	}
	close IN;
	$/="\n";
}

##Note: there are also Pseudo and Undet tRNA genes 
##Note: need to do complement and reverse when caculate codon from anticodon
#############################################
sub parse_tRNA {
	my $file = shift;
	my $data = shift;
	
	open IN,$file || die "fail $file";
	while (<IN>) {
		my @ary = split /\s+/;
		next if($ary[0] eq 'Sequence' || $ary[0] eq 'Name' || $ary[0] eq '--------');
		
		my $seq_id = $ary[0];
		my $tRNA_id = "$ary[0].trna$ary[1]";
		my ($start_pos,$end_pos) = ($ary[2] <= $ary[3]) ? ($ary[2],$ary[3]) : ($ary[3],$ary[2]);
		my $tRNA_strand = ($ary[2] <= $ary[3]) ? '+' : '-';
		my $tRNA_score = $ary[8];
		my $type = $ary[4];
		my $anti_codon = $ary[5];
		
		$data->{$seq_id}{$tRNA_id}{start} = $start_pos;
		$data->{$seq_id}{$tRNA_id}{end} = $end_pos;
		$data->{$seq_id}{$tRNA_id}{score} = $tRNA_score;
		$data->{$seq_id}{$tRNA_id}{strand} = $tRNA_strand;
		$data->{$seq_id}{$tRNA_id}{type} = $type;
		$data->{$seq_id}{$tRNA_id}{anticodon} = $anti_codon;
	}
	close IN;
}


__END__




Input *.tRNA file:

Sequence                tRNA    Bounds  tRNA    Anti    Intron Bounds   Cove
Name            tRNA #  Begin   End     Type    Codon   Begin   End     Score
--------        ------  ----    ------  ----    -----   -----   ----    ------
nscaf1681       1       3208582 3208654 Pseudo  ACT     0       0       21.11
nscaf1681       2       3889539 3889611 Pseudo  ATT     0       0       25.44
nscaf1681       3       4805201 4805275 Pseudo  GTT     0       0       39.59
nscaf1681       4       5387859 5387930 Pseudo  ACT     0       0       25.49
nscaf1681       5       5338374 5338302 Pseudo  GTC     0       0       27.25
nscaf1681       6       2816664 2816592 Pseudo  GTT     0       0       23.20
nscaf1681       7       2238260 2238191 Ser     ACT     0       0       24.44
nscaf1681       8       2134625 2134553 Pseudo  ATC     0       0       20.08
nscaf1681       9       1466253 1466180 Pseudo  ???     0       0       23.66
nscaf1681       10      707463  707390  Pseudo  ATT     0       0       21.87

Input *.tRNA.structure file:

nscaf1071.trna1 (1242603-1242675)       Length: 73 bp
Type: Val       Anticodon: TAC at 34-36 (1242636-1242638)       Score: 76.04
         *    |    *    |    *    |    *    |    *    |    *    |    *    |  
Seq: GGTTCCGTGGTGTAGTGGTtATCACATCTGCTTTACACGCAGAAGGcCGCCAGTTCGATCCTGGCCGGAATCA
Str: >>>>>>>..>>>..........<<<.>>>>>.......<<<<<.....>>>>>.......<<<<<<<<<<<<.

nscaf1071.trna2 (1261594-1261522)       Length: 73 bp
Type: Val       Anticodon: TAC at 34-36 (1261561-1261559)       Score: 76.04
         *    |    *    |    *    |    *    |    *    |    *    |    *    |  
Seq: GGTTCCGTGGTGTAGTGGTtATCACATCTGCTTTACACGCAGAAGGcCGCCAGTTCGATCCTGGCCGGAATCA
Str: >>>>>>>..>>>..........<<<.>>>>>.......<<<<<.....>>>>>.......<<<<<<<<<<<<.

nscaf1108.trna1 (1165590-1165662)       Length: 73 bp
Type: Lys       Anticodon: CTT at 34-36 (1165623-1165625)       Score: 42.49
         *    |    *    |    *    |    *    |    *    |    *    |    *    |  
Seq: CCCCGGCTTGCTCAGTCGGTAGAGCATGAGACTCTTAATCTCAGCGtCGTGGGTTCATGCCCCAGGTCAGGCG
Str: .>>.>>>..>>>>........<<<<.>>>>>.......<<<<<......>>>>.......<<<<.<<<.<<..

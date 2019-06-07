#!/usr/bin/perl

=head1 Program

genemark_convert.pl  --  convert the raw gene predictions to gff, cds, and protein format.

=head1 Description

This program convert the raw result from gene predicting softwares to other 
usually used formats, it is designed to be a universal gene format converter.

The output location files include table format, gff3 format and psl format. The
output sequence files include CDS file, protein file.

There are two options to control gene ID format. The --nametag option will just add a prefix tag 
before the original gene ID; The --finalname option will rename all the genes on the whole genome
level.

There are 4 options to control gene quality. The --perfectgene, --miniseq, --minicds, --checkcds.
In default settings, they are all closed. You can choose one or more of them as your need.


=head1 Command-line Option
  
  $ perl genemark_convert.pl  <predict_result> [genome sequence]
  
  [genome sequence] is used to get cds and protein sequence if <predict_result> doesn't include them.
  --gcode <int>      genetic code, default: 11; supported: 11, 4 and 1
  --nametag <str>    add a prefix tag for gene name on individual sequence level
  --finalname <str>  make final gene name with prefix tag in genome level
  
  --perfectgene      only keep complete genes, remove partial genes
  --miniseq <num>    remove genes predicted on seqeunce shorter than miniseq (default 0 bp)
  --minicds <num>    remove genes which cds length is smaller than minicds (default 0 bp )
  --checkcds         remove genes which is not correct gene model by post-checking on CDS sequence
  --filter_ns        filter CDS which contains non ACGT characters

  --log              transfer STDERR information to a log file
  --verbose          output the running progress information
  --help             output help information

=head1 Usage Exmples

perl predict_convert.pl --final BSP --log --verbose genome.gmhmmp genome.seq

=cut

use strict;
use Getopt::Long;
use Data::Dumper;
use FindBin qw($Bin $Script);

my $Gcode;
my ($Name_tag,$Final_name,$Perfect_gene);
my ($Mini_seq,$Mini_cds,$Check_cds,$Filter_Ns,$Log,$Verbose,$Help);
GetOptions(
	"gcode:i" => \$Gcode,
	"nametag:s" => \$Name_tag,
	"finalname:s"=>\$Final_name,
	"perfectgene" => \$Perfect_gene,
	"miniseq:n" => \$Mini_seq,
	"minicds:n" => \$Mini_cds,
	"checkcds"=>\$Check_cds,
	"filter_ns"=>\$Filter_Ns,
	"log"=>\$Log,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$Gcode ||= 11;
$Name_tag &&= $Name_tag."_";
$Mini_seq ||= 0;
$Mini_cds ||= 0;
die `pod2text $0` if (@ARGV == 0 || $Help);

my $predict_file = shift;
my $sequence_file = shift;


my (%Translate, %hashStarts);
my ($aas, $starts, $base1, $base2, $base3);
if($Gcode == 1)
{
	$aas    = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
	$starts = "---M---------------M---------------M----------------------------";
	$base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
	$base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
	$base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";
}
elsif($Gcode == 4)
{
	$aas    = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
	$starts = "--MM---------------M------------MMMM---------------M------------";
	$base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
	$base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
	$base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";
}
elsif($Gcode == 11)
{
	$aas    = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
	$starts = "---M---------------M------------MMMM---------------M------------";
	$base1  = "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG";
	$base2  = "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG";
	$base3  = "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG";
}
my $temp;
foreach my $i (1..(length $aas))
{
	$temp =substr($base1,($i-1),1);
	$temp.=substr($base2,($i-1),1);
	$temp.=substr($base3,($i-1),1);
	$Translate{$temp}=substr($aas,($i-1),1);
#	if(substr($starts,($i-1),1) eq "M"){
#		$hashStarts{$temp}="M";
#	}
}

my $gff_file = $predict_file.".gff";
my $cds_file = $predict_file.".cds";
my $prot_file = $predict_file.".pep";
my $psl_file = $predict_file.".psl";
my $table_file = $predict_file.".tab";
my $log_file = $predict_file.".log";

$Log ? open(LOG, ">$log_file") : (*LOG = *STDERR);

my %all; ##store the main data
my %seq_len; ##store sequence length
my ($num_multi,$num_single,$num_part,$small_cds_num) = (0,0,0,0); ## store gene number of each type

##reading and parsing the prediction file
read_genemark();
print LOG "\nRead raw result file\n\n" if ($Verbose);

check_redundance(); ## output the checking result to STDERR 

##remove un-wanted genes
perfect_gene() if($Perfect_gene);  ## check the gene type feture

cds_cutoff() if($Mini_cds);

##注意：序列上有问题的基因已经归类于paritial genes
check_cds() if($Check_cds);        ## check cds sequence

filter_Ns() if($Filter_Ns);

##make genome wide gene ID
final_name() if($Final_name);

## output statistic numbers

if ($Verbose && ($Mini_cds || $Check_cds)) {
	print LOG "\n\nStatistic of gene numbers:\n";
	print LOG "\nNote that the genes marked by checkcds are assigned to parital genes \n\n";
	print LOG "  [total genes]   = [multi-exon genes] + [singl-exon genes] + [partial genes]\n";
	print LOG "  [perfect genes] = [multi-exon genes] + [singl-exon genes]\n";
	print LOG "  [final genes]   = [multi-exon genes] + [singl-exon genes] - [small-cds genes]\n\n";
	print LOG "  total genes:      ".($num_multi+$num_single+$num_part)."\n";
	print LOG "  multi-exon genes: $num_multi\n";
	print LOG "  singl-exon genes: $num_single\n";
	print LOG "  partial genes:    $num_part\n";
	print LOG "  perfect genes:    ".($num_multi+$num_single)."\n";
	print LOG "  final genes:      ".($num_multi+$num_single-$small_cds_num)."\n\n";
}

##creating result files
print LOG "Creating result files:\n" if($Verbose);

creat_table();
print LOG "  the result talbe file created\n" if($Verbose);

creat_psl();
print LOG "  the result psl file created\n" if($Verbose);

creat_gff();
print LOG "  the result gff file created\n" if($Verbose);

creat_sequence() if ($sequence_file);
print LOG "  the result sequence files created\n" if($sequence_file && $Verbose);

print LOG "\nWhole task finished\n\n" if($Verbose);


####################################################
################### Sub Routines ###################
####################################################

sub read_genemark{
	my ($seq_name);
	open(PRE,$predict_file) || die("fail to open $predict_file\n");
	while (<PRE>) {
		chomp;
		$seq_name = $1 if (/FASTA definition line: (\S+)\s*/);
		if (/^\s+(\d+)\s+([+-])\s+(\d+)\s+(\d+)/) {
			my (@exon,@orf,@score);
			my ($start,$end) = ($2 eq '+') ? ($3,$4) : ($4,$3);
			push @exon, [$start,$end];
			push @orf, [$start,$end];
			push @score, '.';
			my $gene_name = $seq_name."_orf".make_mark($1,6);
			my $strand = $2;
			my $type = "sigle-exon";
			my $promoter = "none";
			my $polyA = "none";

			$all{$seq_name}{$gene_name}{strand}=$strand;
			$all{$seq_name}{$gene_name}{type}=$type;
			$all{$seq_name}{$gene_name}{promoter}=$promoter;
			$all{$seq_name}{$gene_name}{exon}=\@exon;
			$all{$seq_name}{$gene_name}{orf}=\@orf;
			$all{$seq_name}{$gene_name}{score}=\@score;
			$all{$seq_name}{$gene_name}{polyA}=$polyA;
		}
	}
	close PRE;
	get_seq_len();
}

##from number to string, using as mark
sub make_mark{
	my ($num,$len) = @_;
	my $num_len = length($num);
	if ($num_len < $len) {
		my $zero = '0' x ($len - $num_len);
		$num = $zero.$num;
	}
	return $num;
}


##make genome wide gene ID
sub final_name {
	my $mark = "000001";
	foreach my $seq_name (map {$_->[0]} sort {$a->[1] <=> $b->[1]} map {[$_, /(\d+)$/]} keys %all) {
		my $seq_p = $all{$seq_name};
		foreach my $gene_name (sort keys %$seq_p) {
			my $final_name = $Final_name.$mark;
			$seq_p->{$final_name} = $seq_p->{$gene_name};
			delete $seq_p->{$gene_name};
			$mark++;
		}
	}
}

##creat *.table file
sub creat_table{
	open(OUT,">".$table_file) || die("fail to open $table_file\n");
	foreach my $seq_name (map {$_->[0]} sort {$a->[1] <=> $b->[1]} map {[$_, /(\d+)$/]} keys %all) {
		my $output;
		my $seq_p = $all{$seq_name};
		foreach my $gene_name (sort keys %$seq_p) {
			
			my $gene_p = $seq_p->{$gene_name};
			my $strand = $gene_p->{strand};
			my $type = $gene_p->{type};
			my $promoter = $gene_p->{promoter};
			my ($exon,$exon_num,$cds_size);
			foreach my $p (@{$gene_p->{exon}}) {
				$exon .= "$p->[0]-$p->[1],";
				$exon_num++;
				$cds_size += (abs($p->[1]-$p->[0])+1);
			}
			my $orf;
			foreach my $p (@{$gene_p->{orf}}) {
				$orf .= "$p->[0]-$p->[1],";
			}
			my $score;
			foreach my $value (@{$gene_p->{score}}) {
				$score .= "$value,";
			}
			my $polyA = $gene_p->{polyA};
			$output .= $gene_name."\t".$exon_num."\t".$cds_size."\t".$strand."\t".$seq_name."\t".$seq_len{$seq_name}."\t".$type."\t".$promoter."\t".$polyA."\t".$exon."\t".$orf."\t".$score."\n";
		}

		print OUT $output;
	}
	close OUT;
}

##creat psl file
sub creat_psl{
	my $output;
	#print Dumper \%all;
	foreach my $seq_name (map {$_->[0]} sort {$a->[1] <=> $b->[1]} map {[$_, /(\d+)$/]} keys %all) {
		my $seq_p = $all{$seq_name};
		foreach my $gene_name (sort keys %$seq_p) {
			my $gene_p = $seq_p->{$gene_name};
			my $strand = $gene_p->{strand};
			
			my @oldexon ;
			foreach my $p (@{$gene_p->{exon}}) {
				push @oldexon, [$p->[0],$p->[1]] ;
			}

			if ($strand eq '-') {
				@oldexon = reverse @oldexon;
				foreach my $p (@oldexon) {
					( $p->[0],$p->[1] ) = ( $p->[1],$p->[0] ) ;
				}
			}
			my @exon;
			foreach my $p (@oldexon) {
				push @exon,$p->[0],$p->[1];
			}

			
			my $block = @exon / 2; 
			my ($sizes,$qstarts,$tstarts,$match);
			for (my $i=0; $i<@exon; $i+=2) {
				#print "$exon[$i+1] , $exon[$i]\n";
				$match += $exon[$i+1] - $exon[$i] + 1;
				$sizes .= ( $exon[$i+1] - $exon[$i] + 1 ).",";
				$tstarts .= ( $exon[$i] - 1 ).",";
				$qstarts .= ( $match - ($exon[$i+1] - $exon[$i] + 1) ).",";
			}
			
			my $tstart = $exon[0] - 1;
			my $tend = $exon[(@exon - 1)];

			$output .= "$match\t0\t0\t0\t0\t0\t0\t0\t$strand\t$gene_name\t$match\t0\t$match\t$seq_name\t$seq_len{$seq_name}\t$tstart\t$tend\t$block\t$sizes\t$qstarts\t$tstarts\n";		
			#print "$match\t0\t0\t0\t0\t0\t0\t0\t$strand\t$gene_name\t$match\t0\t$match\t$seq_name\t$seq_len{$seq_name}\t$tstart\t$tend\t$block\t$sizes\t$qstarts\t$tstarts\n";		
		}
	}
	open OUT, ">$psl_file" || die "fail open $psl_file\n";
	print OUT $output;
	close OUT;
}



##creat *.gff file, support gff3
sub creat_gff{
	my $output = "##gff-version 3\n";
	open(OUT,">".$gff_file) || die("fail to open $gff_file\n");
	foreach my $seq_name (map {$_->[0]} sort {$a->[1] <=> $b->[1]} map {[$_, /(\d+)$/]} keys %all) {
		my $seq_p = $all{$seq_name};
		$output .= "##sequence-region $seq_name 1 $seq_len{$seq_name}\n";
		foreach my $gene_name (sort keys %$seq_p) {
			
			my $gene_p = $seq_p->{$gene_name};
			my $strand = $gene_p->{strand};
			#my $type = $gene_p->{type};
			#$my $promoter = $gene_p->{promoter};
			#$my $polyA = $gene_p->{polyA};
			
			##对二维数组操作时要注意：
			## @exon = @{$gene_p->{exon}},只是在一维上新开了@exon,二维仍使用$gene_p->{exon}的内容。
			## 如果不想改变原二维数组内数值顺序，则必须在二维水平上新开@exon数组。
			my (@exon,@orf,@score);
			foreach my $p (@{$gene_p->{exon}}) {
				push @exon,[$p->[0],$p->[1]];
			}
			foreach my $p (@{$gene_p->{orf}}) {
				push @orf,[$p->[0],$p->[1]];
			}
			@score = @{$gene_p->{score}};
			
			my ($gene_start,$gene_end) = ($exon[0][0] < $exon[-1][1]) ? ($exon[0][0], $exon[-1][1]) : ($exon[-1][1], $exon[0][0]);
			
			my $gene_score = (exists $gene_p->{genescore}) ? $gene_p->{genescore} : '.';
			
			my $core_gene_name = $gene_name;
			$core_gene_name =~ s/-\w+$//;
			$output .= "$seq_name\tGeneMark\tgene\t$gene_start\t$gene_end\t.\t$strand\t.\tID=$core_gene_name;Name=$core_gene_name;\n";

			$output .= "$seq_name\tGeneMark\tmRNA\t$gene_start\t$gene_end\t$gene_score\t$strand\t.\tID=$gene_name;Name=$gene_name;Parent=$core_gene_name;\n";
			
			my $mark=@exon;
			$mark=~tr/[0-9]/0/;
			$mark++;
			for (my $i=0; $i<@exon; $i++) {
				my $phase = ($exon[$i][0] - $orf[$i][0]) % 3;
				my ($exon_start,$exon_end) = ($exon[$i][0] < $exon[$i][1]) ? ($exon[$i][0] , $exon[$i][1]) : ($exon[$i][1] , $exon[$i][0]);
				my $cds_id = $gene_name."-C".$mark;
				$output .= "$seq_name\tGeneMark\tCDS\t$exon_start\t$exon_end\t$score[$i]\t$strand\t$phase\tParent=$gene_name;\n";
				$mark++;
			}
		
		}
	}
	print OUT $output;
	close OUT;
}

##check cds on the sequence level, to see frameshift and internal stop codon.
sub check_cds {
	my ($multi_partial_num,$single_partial_num) = (0,0);
	open(IN, $sequence_file) || die ("can not open $sequence_file\n");
	$/=">"; <IN>; $/="\n";
	while (<IN>) {
		my $chr=$1 if(/^(\S+)/);
		$/=">";
		my $seq=<IN>;
		chomp $seq;
		$seq=~s/\s//g;
		$/="\n";

		my $seq_p = $all{$chr};
		foreach my $gene_name (sort keys %$seq_p) {
			my $gene_p = $seq_p->{$gene_name};
			my $strand = $gene_p->{strand};
			my $type = $gene_p->{type};

			my @exon = @{$gene_p->{exon}};
			my $cds_str;
			foreach my $p (@exon) {
				my $str_len = abs($p->[0] - $p->[1]) + 1;
				my $str_start = ($p->[0] < $p->[1]) ?  $p->[0] : $p->[1];
				my $str = substr($seq,$str_start-1,$str_len);
				$str = Complement_Reverse($str) if($strand eq '-');
				$cds_str .= $str;
			}
			
			if( ! check_CDS($cds_str) ){
				$multi_partial_num++ if($type eq "multi-exon");
				$single_partial_num++ if($type eq "sigle-exon");
				$gene_p->{type} = "no-";
				delete $seq_p->{$gene_name};
			}
		}
	}
	print LOG "check cds model on the sequence:\n";
	print LOG "  remove wrong multi-exon genes $multi_partial_num\n";
	print LOG "  remove wrong single-exon genes $single_partial_num\n\n";

	$num_multi -= $multi_partial_num;
	$num_single -= $single_partial_num;
	$num_part += $multi_partial_num + $single_partial_num;
}


##filter Ns in CDS, default 10Ns
sub filter_Ns {
	open(IN, $sequence_file) || die ("can not open $sequence_file\n");
	$/=">"; <IN>; $/="\n";
	while (<IN>) {
		my $chr=$1 if(/^(\S+)/);
		$/=">";
		my $seq=<IN>;
		chomp $seq;
		$seq=~s/\s//g;
		$/="\n";

		my $seq_p = $all{$chr};
		foreach my $gene_name (sort keys %$seq_p) {
			my $gene_p = $seq_p->{$gene_name};
			my $strand = $gene_p->{strand};
			my $type = $gene_p->{type};

			my @exon = @{$gene_p->{exon}};
			my $cds_str;
			foreach my $p (@exon) {
				my $str_len = abs($p->[0] - $p->[1]) + 1;
				my $str_start = ($p->[0] < $p->[1]) ?  $p->[0] : $p->[1];
				my $str = substr($seq,$str_start-1,$str_len);
				$str = Complement_Reverse($str) if($strand eq '-');
				$cds_str .= $str;
			}
			
			my $N_num = $cds_str =~ tr/NnXx//;
			if($N_num > 10){
				delete $seq_p->{$gene_name};
			}
		}
	}
}

#check whether a sequence accord with gene model
#############################################
sub check_CDS{
	my $seq=shift;
	my ($start,$end,$mid,$triple);
	$mid=1;
	my $len=length($seq);
	$triple=1 if($len%3 == 0);
	$start=1 if($seq=~/^ATG/);
	$end=1 if($seq=~/TAA$/ || $seq=~/TAG$/ || $seq=~/TGA$/);
	for (my $i=3; $i<$len-3; $i+=3) {
		my $codon=substr($seq,$i,3);
		$mid=0 if($codon eq 'TGA' || $codon eq 'TAG' || $codon eq 'TAA');
	}
	if ($start && $mid && $end && $triple ) {
		return 1;
	}else{
		return 0;
	}
}


##creat cds and pep file, also exon dna and exon pep file
sub creat_sequence{
	open OUT, ">$cds_file" || die "fail open $cds_file\n";
	open OUT2, ">$prot_file" || die "fail open $prot_file\n";

	open(IN, $sequence_file) || die ("can not open $sequence_file\n");
	$/=">"; <IN>; $/="\n";
	while (<IN>) {
		my $chr=$1 if(/^(\S+)/);
		$/=">";
		my $seq=<IN>;
		chomp $seq;
		$seq=~s/\s//g;
		$/="\n";
		
		my $output;
		my $output2;
		my $output3;
		my $output4;
		my $output5;

		my $seq_p = $all{$chr};
		foreach my $gene_name (sort keys %$seq_p) {
			my $gene_p = $seq_p->{$gene_name};
			my $strand = $gene_p->{strand};
			my $genetype = $gene_p->{type};
			my $gene_name_core = $gene_name;
			$gene_name_core =~ s/-T\w$//; ##real gene name
			my $prot_name = $gene_name;
			$prot_name =~ s/-T(\w)$/-P\1/; ##protein name

			##make gene cds and protein
			my @exon = @{$gene_p->{exon}};
			my @orf = @{$gene_p->{orf}};
			my @score = @{$gene_p->{score}};
			my $exon_num = @exon;
			my ($gene_start,$gene_end)=($exon[0][0] < $exon[-1][1]) ? ($exon[0][0], $exon[-1][1]) : ($exon[-1][1], $exon[0][0]);
			my $cds_str;
			foreach my $p (@exon) {
				
				##染色体从复制起点开始，不会破坏正常的基因结构，因而其实没必要采用环状形式
				##deal with glimmer/genemark circle DNA molecular
				if ($strand eq '+' && $exon[0][0] > $exon[0][1]) {
					my $first_half = substr($seq,$exon[0][0] - 1);
					my $second_half = substr($seq,0,$exon[0][1]);
					$cds_str = $first_half.$second_half;
					last;
				}
				##deal with glimmer/genemark circle DNA molecular
				if ($strand eq '-' && $exon[0][0] < $exon[0][1]) {
					my $first_half = substr($seq,$exon[0][1] - 1);
					my $second_half = substr($seq,0,$exon[0][0]);
					$cds_str = $first_half.$second_half;
					$cds_str = Complement_Reverse($cds_str);
					last;
				}
				
				my $str_len = abs($p->[0] - $p->[1]) + 1;
				my $str_start = ($p->[0] < $p->[1]) ?  $p->[0] : $p->[1];
				my $str = substr($seq,$str_start-1,$str_len);
				$str = Complement_Reverse($str) if($strand eq '-');
				$cds_str .= $str;
			}
			
			my $disp_str;
			Disp_seq(\$cds_str,\$disp_str);
			
			$output .= ">$gene_name  locus=$chr:$gene_start:$gene_end:$strand\n$disp_str";
			
			my $cds_phase = ($orf[0][0] ne 'none') ?  (abs($exon[0][0] - $orf[0][0]) % 3) : 0;
			my $prot_str = cds2aa($cds_str,$cds_phase);
			my $disp_str;
			Disp_seq(\$prot_str,\$disp_str);
			
			$output2 .= ">$prot_name  locus=$chr:$gene_start:$gene_end:$strand\n$disp_str";

			##make exon dna and protein
			my $mark=@exon;
			$mark=~tr/[0-9]/0/;
			$mark++;
			for(my $i=0; $i<@exon; $i++) {
				my $p = $exon[$i];
				my $str;
				my $str_len = abs($p->[0] - $p->[1]) + 1;
				my ($str_start,$str_end) = ($p->[0] < $p->[1]) ?  ($p->[0] , $p->[1]) : ($p->[1] , $p->[0]);
				$str = substr($seq,$str_start-1,$str_len);
				
				$str = Complement_Reverse($str) if($strand eq '-');
				my $disp_str;
				Disp_seq(\$str,\$disp_str);
				my $exon_dna_name = $gene_name.'-E'.$mark;
				my $exon_prot_name = $prot_name.'-E'.$mark;
				$output3 .= ">$exon_dna_name  locus=$chr:$str_start:$str_end:$strand\n$disp_str";
				
				my $exon_phase =($orf[$i][0] ne 'none') ?  (abs($exon[$i][0]-$orf[$i][0]) % 3) : 0;
				my $exon_prot = cds2aa($str,$exon_phase);
				my $disp_str;
				Disp_seq(\$exon_prot,\$disp_str);
				$output4 .= ">$exon_prot_name  locus=$chr:$str_start:$str_end:$strand\n$disp_str";
				
				$mark++;
			}

			##make intron dna file
			my $mark=@exon;
			$mark=~tr/[0-9]/0/;
			$mark++;
			for(my $i=0; $i<@exon-1; $i++) {

				my ($str_start,$str_end) = ($strand eq '+') ?  ($exon[$i][1]+1,$exon[$i+1][0]-1) : ($exon[$i+1][0]+1,$exon[$i][1]-1);
				
				my $str = substr($seq,$str_start-1,$str_end-$str_start+1);
				
				$str = Complement_Reverse($str) if($strand eq '-');
				my $disp_str;
				Disp_seq(\$str,\$disp_str);
				my $intron_dna_name = $gene_name.'-I'.$mark;
				$output5 .= ">$intron_dna_name  locus=$chr:$str_start:$str_end:$strand\n$disp_str";
				
				$mark++;
			}
		}
		
		print OUT $output;
		print OUT2 $output2;
	}

	close(IN);
	close OUT;
	close OUT2;
}

##check whether there is redundance left
##use globle virables
sub check_redundance{
	my ($all_size,$pure_size,$redunt_size,$total_gene_len);
	foreach my $seq_name (sort keys %all) {
		my $seq_p = $all{$seq_name};
		my (@pos1,@pos2);
		foreach my $gene_name (sort keys %$seq_p) {
			my $gene_p = $seq_p->{$gene_name};
			my $strand = $gene_p->{strand};
			my $gene_start = $gene_p->{exon}[0][0];
			my $gene_end = $gene_p->{exon}[-1][1];
			my $gene_len = abs($gene_end-$gene_start) + 1;
			$total_gene_len += $gene_len; ##note
			push @pos1, [$gene_start,$gene_end] if($strand eq '+');
			push @pos2, [$gene_end,$gene_start] if($strand eq '-');
		}
		
		my ($num1,$num2,$num3) = Conjoin_fragment(\@pos1);
		$all_size += $num1;
		$pure_size += $num2;
		$redunt_size += $num3;

		my ($num1,$num2,$num3) = Conjoin_fragment(\@pos2);
		$all_size += $num1;
		$pure_size += $num2;
		$redunt_size += $num3;
	}
	print LOG "\ncheck redundance on bp level:\n  $all_size (all) = $pure_size (pure) + $redunt_size (redunt)\n\n" if($Verbose);
}


##perfect gene set by removing imcomplete genes
##use globle virables
sub perfect_gene{
	my $partial_num;
	foreach my $seq_name (sort keys %all) {
		my $seq_p = $all{$seq_name};
		foreach my $gene_name (sort keys %$seq_p) {
			my $gene_p = $seq_p->{$gene_name};
			my $type = $gene_p->{type};
			if($type =~ /^no-/){
				delete $seq_p->{$gene_name};
				$partial_num++;
			}
		}
	}
	print LOG "perfecting the gene set:\n  remove $partial_num incomplete genes\n\n" if($Verbose);
}

##perfect gene set by removing small genes
##use globle virables
sub cds_cutoff{
	foreach my $seq_name (sort keys %all) {
		my $seq_p = $all{$seq_name};
		foreach my $gene_name (sort keys %$seq_p) {
			my $gene_p = $seq_p->{$gene_name};
			my $cds_size;
			foreach my $p (@{$gene_p->{exon}}) {
				$cds_size += abs($p->[1] - $p->[0]) + 1;
			}
			if($cds_size < $Mini_cds){
				delete $seq_p->{$gene_name};
				$small_cds_num++;
			}
		}
	}
	print LOG "Remove cds less than $Mini_cds bp\n remove $small_cds_num small cds genes\n\n" if($Verbose);
}

sub get_seq_len {
	##get original sequence length
	%seq_len = ();
	open(IN, $sequence_file) || die ("original sequence as second input file is needed for option --backcoor\n");
	$/=">"; <IN>; $/="\n";
	while (<IN>) {
		my $chr=$1 if(/^(\S+)/);
		$/=">";
		my $seq=<IN>;
		chomp $seq;
		$seq=~s/\s//g;
		$/="\n";
		$seq_len{$chr}=length($seq);
	}
	close(IN);
}

#usage: disp_seq(\$string,$num_line);
#############################################
sub Disp_seq{
	my $seq_pp=shift;
	my $disp_pp=shift;
	my $num_line=(@_) ? shift : 50;
	
	my $len=length($$seq_pp);
	for (my $i=0; $i<$len; $i+=$num_line) {
		my $sub=substr($$seq_pp,$i,$num_line);
		$$disp_pp .= $sub."\n";
	}
	$$disp_pp = "\n" if(! $$disp_pp);
}
#############################################


#############################################
sub Complement_Reverse{
	my $seq=shift;
	$seq=~tr/AGCTagct/TCGAtcga/;
	$seq=reverse($seq);
	return $seq;
}
#############################################


## convert cds to protein
sub cds2aa{
	my $seq = shift;
	my $phase = shift || 0;
	die "phase is $phase, it should be 0 1 or 2" if($phase != 0 && $phase != 1 && $phase != 2);

	$seq =~ s/\s//g;
	$seq = uc($seq);
	
	my $len = length($seq);

	my $prot;	
	my $codon = substr($seq,$phase,3);
	for (my $i=$phase; $i<$len; $i+=3) {
		$codon = substr($seq,$i,3);
		last if(length($codon) < 3);
		$prot .= (exists $Translate{$codon}) ? $Translate{$codon} : 'X';
	}
	$prot =~ s/\*$//;
	return $prot;
}


##conjoin the overlapped fragments, and caculate the redundant size
##usage: conjoin_fragment(\@pos);
##		 my ($all_size,$pure_size,$redunt_size) = conjoin_fragment(\@pos);
sub Conjoin_fragment{
	my $pos_p = shift; ##point to the two dimension input array
	my $new_p = [];         ##point to the two demension result array
	
	my ($all_size, $pure_size, $redunt_size) = (0,0,0); 
	
	return (0,0,0) unless(@$pos_p);

	foreach my $p (@$pos_p) {
			($p->[0],$p->[1]) = ($p->[0] <= $p->[1]) ? ($p->[0],$p->[1]) : ($p->[1],$p->[0]);
			$all_size += abs($p->[0] - $p->[1]) + 1;
	}
	
	@$pos_p = sort {$a->[0] <=>$b->[0]} @$pos_p;
	push @$new_p, (shift @$pos_p);
	
	foreach my $p (@$pos_p) {
			if ( ($p->[0] - $new_p->[-1][1]) <= 0 ) { # conjoin
					if ($new_p->[-1][1] < $p->[1]) {
							$new_p->[-1][1] = $p->[1]; 
					}
			}else{  ## not conjoin
					push @$new_p, $p;
			}
	}
	@$pos_p = @$new_p;

	foreach my $p (@$pos_p) {
			$pure_size += abs($p->[0] - $p->[1]) + 1;
	}
	
	$redunt_size = $all_size - $pure_size;
	return ($all_size,$pure_size,$redunt_size);
}

__END__

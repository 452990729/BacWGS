#! /usr/bin/perl

## v2 was modified by liangshuqing at Fri Mar 21 11:39:20 CST 2014
##aim to set the precision of the statistics result
##after running this program, you can get the statistics information on the sreen and the files in separate folders.
## v3: added script used by system into this csript by lss at 201603

=head1
	Usage:
		perl stat.pl [options] genome.seq
	Options:
		--SpeType <B/F/V/P>             set the species type, B: Bacteria, F: Fungi, V: virus, P: phage. [B]	
		--prefix <str>                  set the keyname of result files, needed
		--stat <file>                   set the total stat output file name
		--gene <file>                   set gene cds file in fa format
		--repeat <str>                  stat repeat dir
			--repbase <file>            set repbase result in gff3 format
			--trf <file>                set trf result in gff3 format
		--ncRNA <str>                   stat ncRNA dir
			--rRNA_denovo <file>        set denovo rRNA gff file
			--rRNA_homology <file>      set homology rRNA tab file
			--tRNA <file>               set tRNA gff3 file
			--sRNA <file>               set sRNA gff3 file
			--snRNA <file>              set snRNA gff3 file
			--miRNA <file>              set miRNA gff3 file
		--help                          get the help information

	Version: 2.0
	Contact: liangshuqing@novogene.com

=cut


use strict;
use warnings;

use Getopt::Long;
use File::Basename qw(basename dirname);
use FindBin qw($Bin);
use lib "$Bin/";
use PGAP qw(comma_add);
use File::Path;  ## function " mkpath" and "rmtree" deal with directory

my $SpeType = "B";
my ($Stat,$Gene,$Repeat,$RepBase,$TRF);
my ($NcRNA,$TRNA,$RRNA_denovo,$RRNA_homology,$SRNA,$SnRNA,$MiRNA);
my ($Prefix,$Help);

GetOptions(
	"SpeType:s"=>\$SpeType,
	"stat:s"=>\$Stat,
	"gene:s"=>\$Gene,
	"repeat:s"=>\$Repeat,
	"repbase:s"=>\$RepBase,
	"trf:s"=>\$TRF,
	"ncRNA:s"=>\$NcRNA,
	"tRNA:s"=>\$TRNA,
	"rRNA_denovo:s"=>\$RRNA_denovo,
	"rRNA_homology:s"=>\$RRNA_homology,
	"sRNA:s"=>\$SRNA,
	"snRNA:s"=>\$SnRNA,
	"miRNA:s"=>\$MiRNA,
	"prefix:s"=>\$Prefix,
	"help"=>\$Help
);

die `pod2text $0` if ($Help || @ARGV==0);

#==================Genome===============================
my (%seq);
my (%gene_num,%gene_len);
my (%TRF_num,%repeat_per);
my (%tRNA_num,%rRNA_num,%sRNA_num, %snRNA_num, %miRNA_num);
my ($seq_len, $Genome_gc, $Gene_num, $Gene_len, $Gene_avg_len, $pergenome, $Gene_gc, $Gene_inter_len, $Gene_inter_gc, $Gene_inter_len_per, $total_rRNAd_num, $total_rRNAh_num, $tRNA_sumnum, $sRNA_sumnum, $snRNA_sumnum, $miRNA_sumnum);

$seq_len=0;
##set the statistic precision
my $precision = '%.4f';

my $seq_file = shift;

open SEQ,$seq_file || die "$seq_file $!\n";

$/=">";
while(<SEQ>){
	chomp;
	next if ($_ eq "");
	my($name,$seq)=split(/\n/,$_,2);	
	$name=$1 if ($name=~/^(\S+)/);
	$seq=~s/[\n\r]+//g;
	$seq=uc($seq);
	$seq{$name}=length $seq;
	if (defined $Gene){$gene_num{$name}=0;}else{$gene_num{$name}="--";}
	if (defined $Gene){$gene_len{$name}=0;}else{$gene_len{$name}="--";}
	if (defined $TRF){$TRF_num{$name}=0;}else{$TRF_num{$name}="--";};
	if (defined $TRNA){$tRNA_num{$name}=0;}else{$tRNA_num{$name}="--";}
	if (defined $SRNA){$sRNA_num{$name}=0;}else{$sRNA_num{$name}="--";}
	if (defined $SnRNA){$snRNA_num{$name}=0;}else{$snRNA_num{$name}="--";}
	if (defined $MiRNA){$miRNA_num{$name}=0;}else{$miRNA_num{$name}="--";}
	$seq_len+=length $seq;
}
close SEQ;
$/="\n";

my$gstat=`$Bin/cs $seq_file`;
$Genome_gc=(split(/\n/,$gstat))[-1];
$Genome_gc=(split(/\t/,$Genome_gc))[-1];
$Genome_gc=~s/%//g;
$gstat=~s/\n\t/\n/g;

#===================Gene=============================
if(defined $Gene){
#if(defined $Gene && -d "$All/02.Gene-Prediction/"){
	my$dir=dirname($Gene);
	open OUT,">$dir/$Prefix.gene.stat.xls" || die "$Prefix.gene.stat $!\n";
	open GENE,$Gene || die "$Gene $!\n";
	$Gene_gc=`$Bin/cs $Gene`;
	$Gene_gc=(split(/\n/,$Gene_gc))[-1];
	$Gene_gc=(split(/\t/,$Gene_gc))[-1];
	$Gene_gc=~s/%//g;
	$/=">";
	while(<GENE>){
		chomp;next if ($_ eq "");
		my($head,$seq)=split(/\n/,$_,2);
		my$gene_id=$1 if ($head=~/^(\S+)/);
		my$scaf=$1 if ($head =~ /locus=([\w\d]+)/);
		$seq=~s/[\n\r]+//g;
		$seq=uc($seq);
		$gene_num{$scaf}++;
		$gene_len{$scaf}+=(length $seq);
		$Gene_len += length $seq;
		$Gene_num++;
	}
	$Gene_avg_len = $Gene_len/$Gene_num;
	$pergenome = $Gene_len/$seq_len*100;
	close GENE;
	$Gene_inter_len=$seq_len-$Gene_len;
	$Gene_inter_len_per=$Gene_inter_len*100/$seq_len;
	$Gene_inter_gc=($Genome_gc*$seq_len-$Gene_gc*$Gene_len)/$Gene_inter_len;

	my $temp_seq_len = $seq_len;
	&comma_add (0, $temp_seq_len, $Gene_num, $Gene_len, $Gene_avg_len, $Gene_inter_len);
	&comma_add (2, $Gene_gc, $pergenome, $Gene_inter_gc, $Gene_inter_len_per);
	print OUT "Genome Size:\t$temp_seq_len\n";
	print OUT "Gene Number:\t$Gene_num\n";
	print OUT "Gene Length:\t$Gene_len\n";
	print OUT "GC Content:\t$Gene_gc\n";
	print OUT "\% of Genome(Genes):\t$pergenome\n";
	print OUT "Gene Average Length:\t$Gene_avg_len\n";
	print OUT "Gene Internal Length:\t$Gene_inter_len\n";
	print OUT "Gene Internal GC Content:\t$Gene_inter_gc\n";
	print OUT "% of Genome(internal):\t$Gene_inter_len_per\n";
	close OUT;
	$/="\n";
}

#====================Repeat===========================
my ($TRF_number, $TRF_min, $TRF_max, $TRF_tlen, $TRF_pre);
my ($Min_number, $Min_min, $Min_max, $Min_tlen, $Min_pre);
my ($Mic_number, $Mic_min, $Mic_max, $Mic_tlen, $Mic_pre);
my ($Min_STD_min, $Min_STD_max, $Mic_STD_min, $Mic_STD_max);

#if(defined $Repeat){
if(defined $RepBase) {
        my $Repeat = dirname $RepBase; #added by lss 
		my %hash_RepBase;
		open REPBASE, "$RepBase" || die "$RepBase$!\n";
		while (<REPBASE>) {
			next if (/^\s*#+/);
#Scaffold11 RepeatMasker Transposon 4707 4825 235 - . ID=1_TE01;Target=Sagan-1_AAn 286 398;Class=DNA/TcMar-Sagan;PercDiv=20.1;PercDel=2.5;PercIns=8.0;
			my ($scaf, $sta, $end, $class) = (split (/\t|;/))[0,3,4,10];
			if ($sta > $end) { my $i = $sta; $sta = $end; $end = $i;}
			$class = $1 if ($class =~ /Class=([^\/]+)\//);
			$hash_RepBase{$class}{number} ++;
			$hash_RepBase{Total}{number} ++;
			$hash_RepBase{$class}{len_acc} += $end - $sta + 1;
			$hash_RepBase{Total}{len_acc} += $end - $sta + 1;
			push @{$hash_RepBase{$class}{scaf}{$scaf}}, [$sta, $end];
			push @{$hash_RepBase{Total}{scaf}{$scaf}}, [$sta, $end];
		}
		close REPBASE;
		my $class;
		my %hash_RBout;
		foreach $class (("LTR","DNA","LINE","SINE","RC","Total")) {  #delete "ncRNA" by lss at 20160114
			if (exists $hash_RepBase{$class}) {
				$hash_RBout{$class}{number} = $hash_RepBase{$class}{number};
				$hash_RBout{$class}{avg_len} = $hash_RepBase{$class}{len_acc} / $hash_RepBase{$class}{number};
				foreach my $chr (sort keys %{$hash_RepBase{$class}{scaf}}) {
					$hash_RBout{$class}{total_len} += &Conjoin_fragment($hash_RepBase{$class}{scaf}{$chr});
				}
				$hash_RBout{$class}{pre} = $hash_RBout{$class}{total_len} / $seq_len * 100;
				delete $hash_RepBase{$class};
			}
			else {
				$hash_RBout{$class}{number} = 0;
				$hash_RBout{$class}{avg_len} = 0;
				$hash_RBout{$class}{total_len} = 0;
				$hash_RBout{$class}{pre} = 0;
			}
		}
		if (keys %hash_RepBase > 0) {
			foreach $class (keys %hash_RepBase) {
				$hash_RBout{Unknown}{number} += $hash_RepBase{$class}{number};
				$hash_RBout{Unknown}{len_acc} += $hash_RepBase{$class}{len_acc};
				foreach (keys %{$hash_RepBase{$class}{scaf}}) {
					push @{$hash_RBout{Unknown}{scaf}{$_}}, @{$hash_RepBase{$class}{scaf}{$_}}; 
				}
				delete $hash_RepBase{$class};
			}
			$hash_RBout{Unknown}{avg_len} = $hash_RBout{Unknown}{len_acc} / $hash_RBout{Unknown}{number};
			foreach my $chr (sort keys %{$hash_RBout{Unknown}{scaf}}) {
				$hash_RBout{Unknown}{total_len} += &Conjoin_fragment($hash_RBout{Unknown}{scaf}{$chr});
			}
			$hash_RBout{Unknown}{pre} = $hash_RBout{Unknown}{total_len} / $seq_len * 100;
			delete $hash_RBout{Unknown}{len_acc}; 
			delete $hash_RBout{Unknown}{scaf};
		}
		else {
			($hash_RBout{Unknown}{number}, $hash_RBout{Unknown}{avg_len}, $hash_RBout{Unknown}{total_len}, $hash_RBout{Unknown}{pre}) = (0,0,0,0);
		}
		#output
		open RepBase, ">$Repeat/$Prefix.repbase.stat.xls" || die "can not create file: $Repeat/$Prefix.repbase.stat.xls\n";
		print RepBase "Type\tNumber(#)\tTotal Length(bp)\tIn Genome(%)\tAverage length(bp)\n";
		foreach $class (("LTR","DNA","LINE","SINE","RC","Unknown","Total")) {  #delete ncRNA by lss
			&comma_add (0, $hash_RBout{$class}{number}, $hash_RBout{$class}{avg_len}, $hash_RBout{$class}{total_len});
			&comma_add (4, $hash_RBout{$class}{pre});
			print RepBase "$class\t$hash_RBout{$class}{number}\t$hash_RBout{$class}{total_len}\t$hash_RBout{$class}{pre}\t$hash_RBout{$class}{avg_len}\n";
		}
		close RepBase;
}


if(defined $TRF){
        my $Repeat = dirname($TRF); #added by lss 
		($TRF_number, $TRF_min, $TRF_max, $TRF_tlen, $TRF_pre) = (0, 1000, 0, 0, 0);
		($Min_number, $Min_min, $Min_max, $Min_tlen, $Min_pre) = (0, 1000, 0, 0, 0);
		($Mic_number, $Mic_min, $Mic_max, $Mic_tlen, $Mic_pre) = (0, 1000, 0, 0, 0);
		($Min_STD_min, $Min_STD_max, $Mic_STD_min, $Mic_STD_max) = (10, 60, 2, 6);
		my (%hash_TRF, %hash_Min, %hash_Mic);
		my $out_min_gff = $TRF; $out_min_gff =~ s/trf\.dat\.gff$/Minisatellite.DNA.dat.gff/;
		my $out_mic_gff = $TRF; $out_mic_gff =~ s/trf\.dat\.gff$/Microsatellite.DNA.dat.gff/;
		open MIN, ">$out_min_gff" || die "can not create file: $out_min_gff!\n";
		open MIC, ">$out_mic_gff" || die "can not create file: $out_mic_gff!\n";
		print MIN "##gff-version 3\n";
		print MIC "##gff-version 3\n";
		open TRF,$TRF || die "$TRF $!\n";
		while(my $line = <TRF>){
			chomp $line;
			next if ($line =~ /^#/);
			next if ($line !~ /PeriodSize=(\d+);/);
			my $cur_TRF_len = $1;
			my($scaf,$sta,$end,$copy)=(split(/\s+/, $line))[0,3,4,8];
			if ($sta > $end) { my $i = $sta; $sta = $end; $end = $i;}
			#--all
			$TRF_number ++;
			($cur_TRF_len < $TRF_min) && ($TRF_min = $cur_TRF_len);
			($cur_TRF_len > $TRF_max) && ($TRF_max = $cur_TRF_len);
			push @{$hash_TRF{$scaf}}, [$sta, $end];

			#--Minisatellite DNA
			if ($cur_TRF_len <= $Min_STD_max && $cur_TRF_len >= $Min_STD_min) {
				$Min_number ++;
				print MIN "$line\n";
				($cur_TRF_len < $Min_min) && ($Min_min = $cur_TRF_len);
				($cur_TRF_len > $Min_max) && ($Min_max = $cur_TRF_len);
				push @{$hash_Min{$scaf}}, [$sta, $end];
			}

			#--Microsatellite DNA
			if ($cur_TRF_len <= $Mic_STD_max && $cur_TRF_len >= $Mic_STD_min) {
				$Mic_number ++;
				print MIC "$line\n";
				($cur_TRF_len < $Mic_min) && ($Mic_min = $cur_TRF_len);
				($cur_TRF_len > $Mic_max) && ($Mic_max = $cur_TRF_len);
				push @{$hash_Mic{$scaf}}, [$sta, $end];
			}

			$copy=(split(/;/,$copy))[2];
			$copy=$1 if ($copy =~ /CopyNumber=(\S+)/);
			$TRF_num{$scaf}++;
		}
		close TRF;
		close MIC;
		close MIN;

		$TRF_min = $TRF_max if ($TRF_max == 0);
		$Min_min = $Min_max if ($Min_max == 0);
		$Mic_min = $Mic_max if ($Mic_max == 0);
		foreach my $chr (sort keys %hash_TRF) { $TRF_tlen += &Conjoin_fragment($hash_TRF{$chr});}
		foreach my $chr (sort keys %hash_Min) { $Min_tlen += &Conjoin_fragment($hash_Min{$chr});}
		foreach my $chr (sort keys %hash_Mic) { $Mic_tlen += &Conjoin_fragment($hash_Mic{$chr});}
		$TRF_pre = $TRF_tlen/$seq_len*100;
		$Min_pre = $Min_tlen/$seq_len*100;
		$Mic_pre = $Mic_tlen/$seq_len*100;

		&comma_add (0, $TRF_number, $TRF_min, $TRF_max, $TRF_tlen, $Min_number, $Min_min, $Min_max, $Min_tlen, $Mic_number, $Mic_min, $Mic_max, $Mic_tlen);
		&comma_add (4, $TRF_pre, $Min_pre, $Mic_pre);
		open TRF, ">$Repeat/$Prefix.trf.stat.xls" || die "can not create file: $Repeat/$Prefix.trf.stat.xls!\n";
		print TRF "Type\tNumber(#)\tRepeat Size(bp)\tTotal Length(bp)\tIn Genome(%)\n";
		print TRF "TRF\t$TRF_number\t$TRF_min~$TRF_max\t$TRF_tlen\t$TRF_pre\n";
		print TRF "Minisatellite DNA\t$Min_number\t$Min_min~$Min_max\t$Min_tlen\t$Min_pre\n";
		print TRF "Microsatellite DNA\t$Mic_number\t $Mic_min~$Mic_max\t$Mic_tlen\t$Mic_pre\n";
		close TRF;
}

#===================ncRNA=============================
my($tRNA_len,$rRNA_len,$sRNA_len,$snRNA_len,$miRNA_len)=(0,0,0,0,0);
my($tRNA_per,$rRNA_per,$sRNA_per,$snRNA_per,$miRNA_per)=(0,0,0,0,0);
my($tRNA_avglen,$rRNA_avglen,$sRNA_avglen,$snRNA_avglen,$miRNA_avglen)=(0,0,0,0,0);

if(defined $NcRNA){
	($total_rRNAd_num, $total_rRNAh_num, $tRNA_sumnum, $sRNA_sumnum, $snRNA_sumnum, $miRNA_sumnum) = (0, 0, 0, 0, 0, 0);
	open OUT,">$NcRNA/$Prefix.ncRNA.stat.xls" || die "$!\n";
	print OUT "\tType\tNumber#\tAvg_Len\tTotal_Len\t% in Genome\n";
	if(defined $TRNA){
		read_gff3($TRNA,\$tRNA_sumnum,\$tRNA_avglen,\$tRNA_len,\$tRNA_per,\%tRNA_num);
		&comma_add (0, $tRNA_sumnum, $tRNA_avglen, $tRNA_len);
		&comma_add (4, $tRNA_per);
		print OUT "\ttRNA\t$tRNA_sumnum\t$tRNA_avglen\t$tRNA_len\t$tRNA_per\n";
	}
	else { print OUT "\ttRNA\t-\t-\t-\t-\n"; }
	if(defined $RRNA_denovo){
		my @rRNAd_info;
#		@rRNAd_info = ($SpeType =~ /[BVP]/) ? `perl $Bin/denovo_rRNA_B.pl $RRNA_denovo $seq_file` : `perl $Bin/denovo_rRNA_F.pl $RRNA_denovo $seq_file`;
		@rRNAd_info = ($SpeType =~ /[BVP]/) ? denovo_rRNA_B($RRNA_denovo, $seq_file) : denovo_rRNA_F($RRNA_denovo, $seq_file);
		chomp @rRNAd_info;
		$total_rRNAd_num = pop @rRNAd_info;
		&comma_add (0, $total_rRNAd_num);
		print OUT join ("\n", @rRNAd_info) . "\n";
	}
	else {
		($SpeType =~ /[BVP]/) ? (print OUT "rRNA_de\t5s\t-\t-\t-\t\n", "\t16s\t-\t-\t-\t-\n", "\t23s\t-\t-\t-\t\n") : (print OUT "rRNA_de\t5s\t-\t-\t-\t\n", "\t5.8s\t-\t-\t-\t\n", "\t18s\t-\t-\t-\t-\n", "\t28s\t-\t-\t-\t\n");
	}
	if(defined $RRNA_homology){
		my @rRNAh_info;
#		@rRNAh_info = ($SpeType =~ /[BVP]/) ? `perl $Bin/homology_rRNA_B.pl --rRNA $RRNA_homology --sequence $seq_file` : `perl $Bin/homology_rRNA_F.pl --rRNA $RRNA_homology --sequence $seq_file`;
		@rRNAh_info = ($SpeType =~ /[BVP]/) ? homology_rRNA_B($RRNA_homology,$seq_file) : homology_rRNA_F($RRNA_homology,$seq_file);
		chomp @rRNAh_info;
		$total_rRNAh_num = pop @rRNAh_info;
		&comma_add (0, $total_rRNAh_num);
		print OUT join ("\n", @rRNAh_info) . "\n";
	}
	else {
		($SpeType =~ /[BVP]/) ? (print OUT "rRNA_ho\t5s\t-\t-\t-\t\n", "\t16s\t-\t-\t-\t-\n", "\t23s\t-\t-\t-\t\n") : (print OUT "rRNA_ho\t5s\t-\t-\t-\t\n", "\t5.8s\t-\t-\t-\t\n", "\t18s\t-\t-\t-\t-\n", "\t28s\t-\t-\t-\t\n");
	}
	if(defined $SRNA){
		read_gff3($SRNA,\$sRNA_sumnum,\$sRNA_avglen,\$sRNA_len,\$sRNA_per,\%sRNA_num);
		&comma_add (0, $sRNA_sumnum, $sRNA_avglen, $sRNA_len);
		&comma_add (4, $sRNA_per);
		print OUT "\tsRNA\t$sRNA_sumnum\t$sRNA_avglen\t$sRNA_len\t$sRNA_per\n";
	}
	else { print OUT "\tsRNA\t-\t-\t-\t-\n"; }
	if(defined $SnRNA){
		read_gff3($SnRNA,\$snRNA_sumnum,\$snRNA_avglen,\$snRNA_len,\$snRNA_per,\%snRNA_num);
		&comma_add (0, $snRNA_sumnum, $snRNA_avglen, $snRNA_len);
		&comma_add (4, $snRNA_per);
		print OUT "\tsnRNA\t$snRNA_sumnum\t$snRNA_avglen\t$snRNA_len\t$snRNA_per\n";
	}
	#else { print OUT "\tsnRNA\t-\t-\t-\t-\n"; }
	if(defined $MiRNA){
		read_gff3($MiRNA,\$miRNA_sumnum,\$miRNA_avglen,\$miRNA_len,\$miRNA_per,\%miRNA_num);
		&comma_add (0, $miRNA_sumnum, $miRNA_avglen, $miRNA_len);
		&comma_add (4, $miRNA_per);
		print OUT "\tmiRNA\t$miRNA_sumnum\t$miRNA_avglen\t$miRNA_len\t$miRNA_per\n";
	}
	#else { print OUT "\tmiRNA\t-\t-\t-\t-\n"; }
	close OUT;
}

#========================== stat all result ==========================#
if (defined $Stat) {
	open OUT, ">$Stat" || die;
	&comma_add(0,$seq_len) if (defined $seq_len);
	print OUT "Genome Size(bp):\t", ((defined $seq_len)?$seq_len:"-"), "\n";
	print OUT "GC Content(%):\t", ((defined $Genome_gc)?$Genome_gc:"-"), "\n";
	print OUT "Gene Number(#):\t", ((defined $Gene_num)?$Gene_num:"-"), "\n";
	print OUT "Gene Length(bp):\t", ((defined $Gene_len)?$Gene_len:"-"), "\n";
	print OUT "Gene Average Length(bp):\t", ((defined $Gene_avg_len)?$Gene_avg_len:"-"), "\n";
	print OUT "Gene Length/Genome(%):\t", ((defined $pergenome)?$pergenome:"-"), "\n";
	print OUT "GC Content in Gene Region(%):\t", ((defined $Gene_gc)?$Gene_gc:"-"), "\n";
	print OUT "Intergenic Region Length(bp):\t", ((defined $Gene_inter_len)?$Gene_inter_len:"-"), "\n";
	print OUT "GC Content in Intergenic Region(%):\t", ((defined $Gene_inter_gc)?$Gene_inter_gc:"-"), "\n";
	print OUT "Intergenic Region Length/Genome(%):\t", ((defined $Gene_inter_len_per)?$Gene_inter_len_per:"-"), "\n";
	print OUT "Tandem Repeat Number(#):\t", ((defined $TRF_number)?$TRF_number:"-"), "\n";
	print OUT "Tandem Repeat Length(bp):\t", ((defined $TRF_tlen)?$TRF_tlen:"-"), "\n";
	print OUT "Tandem Repeat Size(bp):\t", ((defined $TRF_min && defined $TRF_max)?"$TRF_min-$TRF_max":"-"), "\n";
	print OUT "Tandem Repeat Length/Genome(%):\t", ((defined $TRF_pre)?$TRF_pre:"-"), "\n";
	print OUT "Minisatellite DNA Number(#):\t", ((defined $Min_number)?$Min_number:"-"), "\n";
	print OUT "Microsatellite DNA Number(#):\t", ((defined $Mic_number)?$Mic_number:"-"), "\n";
	print OUT "rRNA Number(#):\t",((defined $total_rRNAd_num && $total_rRNAd_num ne "0")?$total_rRNAd_num:((defined $total_rRNAh_num)?$total_rRNAh_num:"-")),"\n";
	print OUT "tRNA Number(#):\t", ((defined $tRNA_sumnum)?$tRNA_sumnum:"-"), "\n";
	print OUT "sRNA Number(#):\t", ((defined $sRNA_sumnum)?$sRNA_sumnum:"-"), "\n";
	print OUT "snRNA Number(#):\t", ((defined $snRNA_sumnum)?$snRNA_sumnum:"-"), "\n";
	print OUT "miRNA Number(#):\t", ((defined $miRNA_sumnum)?$miRNA_sumnum:"-"), "\n";
	print OUT "Genomic Island Number(#):\t-\n";
	print OUT "Prophage Number(#):\t-\n";
	close OUT;
}


#==================sub program========================
sub read_gff3{
	my($gff3,$sum_num,$avg_len,$len,$percentage,$num)=@_;
	die "CHECK: $gff3 is not exists !" if (!-e $gff3);
	open GFF,$gff3 || die "$gff3 $!";
	while(<GFF>){
		next if (/^#/);
		my$scaf=(split)[0];
		$$sum_num++;
		$num->{$scaf}++;
	}
#	$$len=`perl  $Bin/stat_TE.pl -gff $gff3 -rank all`;
	$$len= stat_TE("$gff3","all");

	$$len=$1 if ($$len=~/(\d+)/);
	$$len=0 if ($$len eq "");
	if($$sum_num != 0){$$avg_len = sprintf("%d", $$len/$$sum_num);}else{$$avg_len=0;}
	$$percentage = $$len/$seq_len*100;
	close GFF;
}
sub Conjoin_fragment{
	my $pos_p = shift; ##point to the two dimension input array
	my $distance = shift || 0;
	my $new_p = [];         ##point to the two demension result array

	my ($all_size, $pure_size, $redunt_size) = (0,0,0);

	return (0,0,0) unless(@$pos_p);

	foreach my $p (@$pos_p) {
		($p->[0],$p->[1]) = ($p->[0] <= $p->[1]) ? ($p->[0],$p->[1]) : ($p->[1],$p->[0]);
		$all_size += abs($p->[0] - $p->[1]) + 1;
	}

	@$pos_p = sort {$a->[0] <=> $b->[0]} @$pos_p;
	push @$new_p, (shift @$pos_p);

	foreach my $p (@$pos_p) {
		if ( ($p->[0] - $new_p->[-1][1]) <= $distance ) { # conjoin two neigbor fragements when their distance lower than 10bp
			if ($new_p->[-1][1] < $p->[1]) { $new_p->[-1][1] = $p->[1];}
		}else{ push @$new_p, $p;} ## not conjoin
	}
	@$pos_p = @$new_p;

	foreach my $p (@$pos_p) { $pure_size += abs($p->[0] - $p->[1]) + 1;}

	$redunt_size = $all_size - $pure_size;
	return $pure_size;
}
#===================================from origin script: stat_TE.pl
sub stat_TE{
	my ($file,$Rank) = @_;
	my %Data;
#	foreach my $file (@Gff) {
	Read_gff($file,\%Data,$Rank);
#	}
	foreach my $TE_type (sort keys %Data) {
		my $TE_type_p = $Data{$TE_type};
		my $TE_type_size = 0;
		foreach my $chr (sort keys %$TE_type_p) {
			$TE_type_size += Conjoin_fragment($TE_type_p->{$chr});
		}
		return "$TE_type\t$TE_type_size\n";
	}
}
sub Read_gff{
	my $file=shift;
	my $hash_p=shift; 
	my $Rank=shift;

	open (IN,$file) || die ("fail open $file\n");
	while (<IN>) {
		chomp;
		s/^\s+//;
		next if(/^\#/);
		my @temp=split(/\t/);

		my $tname = $temp[0];
		my $strand = $temp[6];
		my $start = $temp[3];
		my $end = $temp[4];
		my $TE_name = $1 if ($temp[8] =~ /Target=([^; ]+);*/);
		my $TE_class = $1 if ($temp[8] =~ /Class=([^;]+);*/);
		$TE_class && $TE_class =~ s/\?$//;
		if ($TE_class) {next if($TE_class eq "scRNA" || $TE_class eq "srpRNA")};
		if ($TE_class) {$TE_class = "Other/DNA_virus" if($TE_class eq "DNA_virus")};


		my ($TE_type,$TE_subtype);
		if ($TE_class) {$TE_type = $1 if($TE_class =~ /^([^\/]+)\/*/)};
		if ($TE_class) {$TE_type = "DNA" if($TE_type eq "RC")};
		if ($TE_class) {$TE_subtype = $1 if($TE_class =~ /^[^\/]+\/([^\/]+)/)};
		$TE_subtype ||= $TE_type;

		my $TE_supertype = $TE_subtype; ##±£»¤ÒÔÏÂ¼¸Àà
#		if($TE_subtype ne "En-Spm" && $TE_subtype ne "Dong-R4" && $TE_subtype ne "Rex-Babar" && $TE_subtype ne "Y-chromosome"){
#			$TE_supertype =~ s/-[^-]+$//; ##É¾³ý-ºóÃæ×ÓÀà²¿·Ö
#			$TE_supertype =~ s/\(.+\)$//; ##É¾³ý´øÓÐÀ¨ºÅµÄ×ÓÀà²¿·Ö
#		} by lss at 201603

		my $need_stat;
		if ($Rank eq "all") {
			$need_stat = "Total_TE";
		}elsif($Rank eq "type"){
			$need_stat = $TE_type;
		}elsif($Rank eq "supertype"){
			$need_stat = "$TE_type/$TE_supertype";
		}elsif($Rank eq "subtype"){
			$need_stat = "$TE_type/$TE_subtype";
		}elsif($Rank eq "family"){
			$need_stat = "$TE_type/$TE_subtype/$TE_name";
		}

		push @{$hash_p->{$need_stat}{$tname}}, [$start,$end]; 
	}
	close(IN);
}
#===================================END of script: stat_TE.pl
#===================================From script: denovo_rRNA_B.pl
#@rRNAd_info = ($SpeType =~ /[BVP]/) ? denovo_rRNA_B($RRNA_denovo $seq_file) : denovo_rRNA_F($RRNA_denovo $seq_file);
sub denovo_rRNA_B{
	my ($rRNA, $sequence) = @_;
	my $rRNA_num = 0;
	open RNA,$rRNA or die 'cannot open the file, ', "$!\n";
	my ($count_5s, $count_16s, $count_23s);
	my ($length_5s, $length_16s, $length_23s);
	while(<RNA>){
		next if(/^#/);
		my ($seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attribute) = split;
		if($attribute eq '5s_rRNA'){
			++$count_5s;
			$length_5s += ($end-$start+1);
		}
		elsif($attribute eq '16s_rRNA'){
			++$count_16s;
			$length_16s += ($end-$start+1);
		}elsif($attribute eq '23s_rRNA'){
			++$count_23s;
			$length_23s += ($end-$start+1);
		}else{
			warn("unknown attribute! $attribute in file: $rRNA\n");
		}	
	}
	close RNA;

	my $total_length = &readfasta($sequence);
	$length_5s ||= 0;
	$length_16s ||= 0;
	$length_23s ||= 0;
	my $total_ratio = sprintf("%.4f", 100*($length_5s+$length_16s+$length_23s)/$total_length);
#print "Type\tCopy#\tAvg_Len\tTotal_Len\t% in Genome\n";
#print "rRNA_de";
	my @out;
	if (defined $count_5s){
		my $temp_avg = $length_5s/$count_5s;
		&comma_add (0, $count_5s, $temp_avg, $length_5s);
		push @out, "rRNA_de\t5s\t$count_5s\t$temp_avg\t$length_5s\n";
		$rRNA_num += $count_5s;
	}else{
		push @out,"rRNA_de\t5s\t0\t0\t0";
	}
	if (defined $count_16s){
		my $temp_avg = $length_16s/$count_16s;
		&comma_add (0, $count_16s, $temp_avg, $length_16s);
		&comma_add (4, $total_ratio);
		push @out, "\t16s\t$count_16s\t$temp_avg\t$length_16s\t$total_ratio\n";
		$rRNA_num += $count_16s;
	}else{
		&comma_add (4, $total_ratio);
		push @out, "\t16s\t0\t0\t0\t$total_ratio";
	}
	if (defined $count_23s){
		my $temp_avg = $length_23s/$count_23s;
		&comma_add (0, $count_23s, $temp_avg, $length_23s);
		push @out,"\t23s\t$count_23s\t$temp_avg\t$length_23s\t";
		$rRNA_num += $count_23s;
	}else{
		push @out, "\t23s\t0\t0\t0";
	}

	push @out, "$rRNA_num";
	return @out;
}

sub readfasta{
	my $fa_file = shift or die 'no parameter input ', "\n";
	my $fa_length = 0;
	my $original = $/;
	open FASTA, $fa_file or die 'cannot open the file, ', $!, ". \n";
	$/ = '>';
	while(<FASTA>){
		chomp;
		next unless $_;
		my ($id, $seq) = split(/\n/, $_, 2);
		$seq =~ s/[\n\s]//mg;
		$fa_length += length($seq);
	}
	$/ = $original;
	close FASTA;
	return $fa_length;
}

#sub comma_add {
#	my $nu = shift;
#	my $arg = "%.${nu}f";
#	foreach (@_) {
#		$_ = sprintf($arg,$_);
#		$_ = /(\d+)\.(\d+)/ ? comma($1) . '.' . $2 : comma($_);
#		$_ = "0" if ($_ =~ /^[0\.]+$/);
#	}
#}
sub comma{
	my ($c,$rev) = @_;
	(length($c) > 3) || return($c);
	$rev || ($c = reverse $c);
	$c =~ s/(...)/$1,/g;
	$rev || ($c = reverse $c);
	$rev ? ($c =~ s/,$//) : ($c =~ s/^,//);
	$c;
}
#===================================END of script: denovo_rRNA_B.pl
#===================================from script: denovo_rRNA_F.pl
sub denovo_rRNA_F{
	my ($rRNA, $sequence) = @_;
	my ($count_5s, $count_58s, $count_18s, $count_28s);
	my ($length_5s, $length_58s, $length_18s, $length_28s);

	my $rRNA_num = 0;

	open RNA,$rRNA or die 'cannot open the file, ', "$!\n";
	while(<RNA>){
		next if(/^#/);
		my ($seqname, $source, $feature, $start, $end, $score, $strand, $frame, $attribute) = split;
		if($attribute eq '5s_rRNA')	{
			++$count_5s;
			$length_5s += ($end-$start+1);
		}elsif($attribute eq '5.8s_rRNA')	{
			++$count_58s;
			$length_58s += ($end-$start+1);
		}elsif($attribute eq '18s_rRNA')	{
			++$count_18s;
			$length_18s += ($end-$start+1);
		}elsif($attribute eq '28s_rRNA')	{
			++$count_28s;
			$length_28s += ($end-$start+1);
		}else	{
			warn("unknown attribute! $attribute in file: $rRNA\n");
		}	
	}
	close RNA;

	my $total_length = &readfasta($sequence);
	$length_5s ||= 0;
	$length_58s ||= 0;
	$length_18s ||= 0;
	$length_28s ||= 0;
	my $total_ratio = sprintf("%.4f", 100*($length_5s+$length_58s+$length_18s+$length_28s)/$total_length);
	&comma_add (4, $total_ratio);
	#print "Type\tCopy#\tAvg_Len\tTotal_Len\t% in Genome\n";
	#print "rRNA_de";
	my @out;
	if (defined $count_5s){
		my $temp_avg = $length_5s/$count_5s;
		&comma_add (0, $count_5s, $temp_avg, $length_5s);
		push @out, "rRNA_de\t5s\t$count_5s\t$temp_avg\t$length_5s";
		$rRNA_num += $count_5s;
	}else{
		push @out, "rRNA_de\t5s\t0\t0\t0";
	}
	if (defined $count_58s){
		my $temp_avg = $length_58s/$count_58s;
		&comma_add (0, $count_58s, $temp_avg, $length_58s);
		push @out, "\t5.8s\t$count_58s\t$temp_avg\t$length_58s";
		$rRNA_num += $count_58s;
	}else{
		push @out, "\t5.8s\t0\t0\t0";
	}
	if (defined $count_18s){
		my $temp_avg = $length_18s/$count_18s;
		&comma_add (0, $count_18s, $temp_avg, $length_18s);
		push @out, "\t18s\t$count_18s\t$temp_avg\t$length_18s\t$total_ratio";
		$rRNA_num += $count_18s;
	}else{
		push @out, "\t18s\t0\t0\t0\t$total_ratio"
	}
	if (defined $count_28s){
		my $temp_avg = $length_28s/$count_28s;
		&comma_add (0, $count_28s, $temp_avg, $length_28s);
		push @out, "\t28s\t$count_28s\t$temp_avg\t$length_28s";
		$rRNA_num += $count_28s;
	}else{
		push @out, "\t28s\t0\t0\t0";
	}

	push @out, "$rRNA_num";
	return @out;
}
#===================================END of script: denovo_rRNA_F.pl
#===================================From script: homology_rRNA_B.pl
##@rRNAh_info = ($SpeType =~ /[BVP]/) ? `perl $Bin/homology_rRNA_B.pl --rRNA $RRNA_homology --sequence $seq_file` : `perl $Bin/homology_rRNA_F.pl --rRNA $RRNA_homology --sequence $seq_file`;
sub homology_rRNA_B{
	my ($file_path, $sequence) = @_;
	my ($count_5s, $count_16s, $count_23s);
	my ($length_5s, $length_16s, $length_23s);
	my %rRNA;
	my $rRNA_num = 0;
	my $cutoff = 0.8;

	open TABLE, $file_path or die 'cannot open the gff file', "$!\n";
	while(<TABLE>)	{
		next if(/^#/);
		chomp;
##1:Query_id  2:Query_length  3:Query_start  4:Query_end   5:Subject_id  6:Subject_length  7:Subject_start  8:Subject_end 
##9:Identity  10:Positive  11:Gap  12:Align_length  13:Score  14:E_value  15:Query_annotation  16:Subject_annotation
		my ($query_id, $query_length, $query_start, $query_end, $subject_id, $subject_length, $subject_start, $subject_end, $identity, $positive, $gap, $align_length, $score, $e_value, $query_anno, $subject_anno) = split;
		my $rRNA_type = 'unknown_type';
		next if($identity < 0.01);
		my $threshold = $align_length / $query_length;
		next if ($threshold < $cutoff);
		$positive = ($subject_start > $subject_end)?'-':'+';
		($subject_start, $subject_end) = ($positive eq '+') ? ($subject_start, $subject_end) : ($subject_end, $subject_start);
		if($query_id =~ /\w+\|(\S+)\|\S+\|(\S+)$/)		{
			$rRNA_type = $2;
		}elsif($query_id =~ /\S+\_(\S+)$/)		{
			$rRNA_type = $1;
		}
		next if($rRNA_type eq '-|');
		my @tmp = ($subject_length, $subject_start, $subject_end, $rRNA_type, $positive);
		push(@{$rRNA{$subject_id}}, \@tmp);
	};
	close TABLE;

#print Dumper(\%rRNA);
##if need add output the gff format file, it will be used!
#print "seqname\tsource\tfeature\tstart\tend\tscore\t+/-\tframe\tattribute\n";
	foreach my $scaf (sort {$a cmp $b} (keys %rRNA)){
		my $filter = &merge_rrna($rRNA{$scaf});
		my @filter_rRNA = @{$filter};
		#print Dumper(\@filter_rRNA);
		for(my $index=0; $index<@filter_rRNA; ++$index)	{
			next if (!defined($filter_rRNA[$index]));
			my @tmp = @{$filter_rRNA[$index]};
			$tmp[3] = lc($tmp[3]);
			if($tmp[3] eq '5s'){
				++$count_5s;
				$length_5s += ($tmp[2]-$tmp[1]+1);
			}elsif($tmp[3] eq '16s'){
				++$count_16s;
				$length_16s += ($tmp[2]-$tmp[1]+1);
			}elsif($tmp[3] eq '23s'){
				++$count_23s;
				$length_23s += ($tmp[2]-$tmp[1]+1);
			}else{
				warn("unknown attribute!\n");
			}	
			#print "$scaf\tblastn\trRNA\t$tmp[1]\t$tmp[2]\t$tmp[4]\t--\t--\t$tmp[3]_rRNA\n"
		}
	}
	##
	my $total_length = &readfasta($sequence);
	$length_5s ||= 0;
	$length_16s ||= 0;
	$length_23s ||= 0;
	my $total_ratio = sprintf("%.4f", 100*($length_5s+$length_16s+$length_23s)/$total_length);
	&comma_add (4, $total_ratio);
	##stat the info of rRNA
	##print "Type\tCopy#\tAvg_Len\tTotal_Len\t% in Genome\n";
	#print "rRNA_ho\t";
	my @out;
	if (defined $count_5s){
		my $temp_avg = $length_5s/$count_5s;
		&comma_add (0, $count_5s, $temp_avg, $length_5s);
		push @out,  "rRNA_ho\t5s\t$count_5s\t$temp_avg\t$length_5s";
		$rRNA_num += $count_5s;
	}else{
		push @out,   "rRNA_ho\t5s\t0\t0\t0";
	}
	if (defined $count_16s){
		my $temp_avg = $length_16s/$count_16s;
		&comma_add (0, $count_16s, $temp_avg, $length_16s);
		push @out,  "\t16s\t$count_16s\t$temp_avg\t$length_16s\t$total_ratio";
		$rRNA_num += $count_16s;
	}else{	
		push @out,  "\t16s\t0\t0\t0\t$total_ratio";
	}
	if (defined $count_23s){
		my $temp_avg = $length_23s/$count_23s;
		&comma_add (0, $count_23s, $temp_avg, $length_23s);
		push @out,  "\t23s\t$count_23s\t$temp_avg\t$length_23s";
		$rRNA_num += $count_23s;
	}else{
		push @out,  "\t23s\t0\t0\t0";
	}

	push @out,  "$rRNA_num";
	return @out;
}

sub merge_rrna{
	my $ref = shift;
	my @sub_rRNA = @{$ref};
	my $size = @sub_rRNA;
	foreach (@sub_rRNA)	{
		#	next if(!defined $_);
		#	print "@{$_}", "\n";
	}
	for(my $index=0; $index < $size; ++$index) {
		next if(!defined($sub_rRNA[$index]));
		for(my $loop=0; $loop < $size; ++$loop)	{
			next if ($loop == $index);
			next if (!defined($sub_rRNA[$index]));
			next if (!defined($sub_rRNA[$loop]));
			my @outer = @{$sub_rRNA[$index]};
			my @inner = @{$sub_rRNA[$loop]};
			if($outer[3] eq $inner[3] and $outer[4] eq $inner[4]) {
				if($outer[1] == $inner[1]) {
					$outer[2] >= $inner[2] ? (delete($sub_rRNA[$loop])) : (delete($sub_rRNA[$index]));
				}elsif($outer[1] > $inner[1]){
					if($outer[2] <= $inner[2]){
						delete($sub_rRNA[$index]);
					}else{
						$inner[2] = $outer[2];
					}
				}else{
					if($inner[2] <= $outer[2]){
						delete($sub_rRNA[$loop]);
					}else{
						$outer[2] = $inner[2];
					}
				}
			}
		}
	}
	return \@sub_rRNA;
}
#===================================END of script: homology_rRNA_B.pl
#===================================From script: homology_rRNA_F.pl
sub homology_rRNA_F{
	my ($file_path, $sequence) = @_;#!/usr/bin/perl -w

	my ($count_5s, $count_58s, $count_18s, $count_28s);
	my ($length_5s, $length_58s, $length_18s, $length_28s);
	my %rRNA;
	my $rRNA_num = 0;
	my $cutoff = 0.8;

	open TABLE, $file_path or die 'cannot open the gff file', "$!\n";
	while(<TABLE>){
		next if(/^#/);
		chomp;
##1:Query_id  2:Query_length  3:Query_start  4:Query_end   5:Subject_id  6:Subject_length  7:Subject_start  8:Subject_end 
##9:Identity  10:Positive  11:Gap  12:Align_length  13:Score  14:E_value  15:Query_annotation  16:Subject_annotation
		my ($query_id, $query_length, $query_start, $query_end, $subject_id, $subject_length, $subject_start, $subject_end, $identity, $positive, $gap, $align_length, $score, $e_value, $query_anno, $subject_anno) = split;
		my $rRNA_type = 'unknown_type';
		next if($identity < 0.01);
		my $threshold = $align_length / $query_length;
		next if ($threshold < $cutoff);
		$positive = ($subject_start > $subject_end)?'-':'+';
		($subject_start, $subject_end) = ($positive eq '+') ? ($subject_start, $subject_end) : ($subject_end, $subject_start);
		if($query_id =~ /\w+\|(\S+)\|\S+\|(\S+)$/){
			$rRNA_type = $2;
		}elsif($query_id =~ /\S+\_(\S+)$/){
			$rRNA_type = $1;
		}
		next if($rRNA_type eq '-|');
		my @tmp = ($subject_length, $subject_start, $subject_end, $rRNA_type, $positive);
		push(@{$rRNA{$subject_id}}, \@tmp);
	};
	close TABLE;
	#print Dumper(\%rRNA);
	##if need add output the gff format file, it will be used!
	#print "seqname\tsource\tfeature\tstart\tend\tscore\t+/-\tframe\tattribute\n";
	foreach my $scaf (sort {$a cmp $b} (keys %rRNA)){
		my $filter = &merge_rrna($rRNA{$scaf});
		my @filter_rRNA = @{$filter};
		#print Dumper(\@filter_rRNA);
		for(my $index=0; $index<@filter_rRNA; ++$index){
			next if (!defined($filter_rRNA[$index]));
			my @tmp = @{$filter_rRNA[$index]};
			$tmp[3] = lc($tmp[3]);
			if($tmp[3] eq '5s'){
				++$count_5s;
				$length_5s += ($tmp[2]-$tmp[1]+1);
			}
			if($tmp[3] eq '5.8s'){
				++$count_58s;
				$length_58s += ($tmp[2]-$tmp[1]+1);
			}elsif($tmp[3] eq '18s'){
				++$count_18s;
				$length_18s += ($tmp[2]-$tmp[1]+1);
			}elsif($tmp[3] eq '28s'){
				++$count_28s;
				$length_28s += ($tmp[2]-$tmp[1]+1);
			}else{
				warn("unknown attribute!\n");
			}	
			#print "$scaf\tblastn\trRNA\t$tmp[1]\t$tmp[2]\t$tmp[4]\t--\t--\t$tmp[3]_rRNA\n"
		}
	}
	##
	my $total_length = &readfasta($sequence);
	$length_5s ||= 0;
	$length_58s ||= 0;
	$length_18s ||= 0;
	$length_28s ||= 0;
	my $total_ratio = sprintf("%.4f", 100*($length_5s+$length_58s+$length_18s+$length_28s)/$total_length);
	&comma_add (4, $total_ratio);
	##stat the info of rRNA
	##print "Type\tCopy#\tAvg_Len\tTotal_Len\t% in Genome\n";
	#print "rRNA_ho\t";
	my @out;
	if (defined $count_5s){
		my $temp_avg = $length_5s/$count_5s;
		&comma_add (0, $count_5s, $temp_avg, $length_5s);
		push @out,  "rRNA_ho\t5s\t$count_5s\t$temp_avg\t$length_5s";
		$rRNA_num += $count_5s;
	}else{
#print  "rRNA_ho\t5s\t0\t0\t0"; # by lss at 201605
		push @out,  "rRNA_ho\t5s\t0\t0\t0";
	}
	if (defined $count_58s){
		my $temp_avg = $length_58s/$count_58s;
		&comma_add (0, $count_58s, $temp_avg, $length_58s);
		push @out,  "\t5.8s\t$count_58s\t$temp_avg\t$length_58s";
		$rRNA_num += $count_58s;
	}else{
		push @out,  "\t5.8s\t0\t0\t0";
	}
	if (defined $count_18s){
		my $temp_avg = $length_18s/$count_18s;
		&comma_add (0, $count_18s, $temp_avg, $length_18s);
		push @out,  "\t18s\t$count_18s\t$temp_avg\t$length_18s\t$total_ratio";
		$rRNA_num += $count_18s;
	}else{
		push @out,  "\t18s\t0\t0\t0\t$total_ratio";
	}
	if (defined $count_28s){
		my $temp_avg = $length_28s/$count_28s;
		&comma_add (0, $count_28s, $temp_avg, $length_28s);
		push @out,  "\t28s\t$count_28s\t$temp_avg\t$length_28s";
		$rRNA_num += $count_28s;
	}else{
		push @out,  "\t28s\t0\t0\t0";
	}
	push @out,  "$rRNA_num";
	return @out;
}
#===================================END of script: homology_rRNA_F.pl

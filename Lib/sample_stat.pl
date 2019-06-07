#!/usr/bin/perl -w
use strict;
(@ARGV >= 3) || die"Usage: perl $0 <ass.list> <01.run_component> <SpeType> [outdir] \n";
my ($ass, $run_dir, $SpeType, $outdir) = @ARGV;
$outdir ||= "./";
(-s $ass) && (-d $run_dir) || die$!;

chomp(my @names =`awk '{print \$1}' $ass`);

##stat for all sample
my ($gene_all,$repbase_all,$trf_all,$ncRNA_all,$transposon_all_nuc,$transposon_all_pep);
foreach my $sample (@names) {
    if (-e "$run_dir/$sample/01.Gene_Prediction/$sample.gene.stat.xls") {
        my @column = `awk -F '\t' '{print \$2}' $run_dir/$sample/01.Gene_Prediction/$sample.gene.stat.xls`;
        map {$_ =~ s/\n//g} @column[0,1,2,5,4];
        $gene_all .=  "$sample\t";
        $gene_all .= join ("\t", @column[0,1,2,5,4]);
        $gene_all .=  "\n";
    }
    if (-e "$run_dir/$sample/02.Repeat-Finding/01.repbase/$sample.repbase.stat.xls") {
        open IN,"$run_dir/$sample/02.Repeat-Finding/01.repbase/$sample.repbase.stat.xls"||die$!;
        <IN>;
        while(<IN>){$repbase_all .= "$sample\t$_"}
        close IN;
    }
    if (-e "$run_dir/$sample/02.Repeat-Finding/02.trf/$sample.trf.stat.xls") {
        open IN,"$run_dir/$sample/02.Repeat-Finding/02.trf/$sample.trf.stat.xls"||die$!;
        <IN>;
        while(<IN>){$trf_all .= "$sample\t$_"}
        close IN;        
    }
   
    if (-e "$run_dir/$sample/03.ncRNA-Finding/$sample.ncRNA.stat.xls") {
        my $i=0; my @sample_ncRNA;
        my ($miRNA, $snRNA) = (0, 0);
        for (`less $run_dir/$sample/03.ncRNA-Finding/$sample.ncRNA.stat.xls`) {
            chomp(@{$sample_ncRNA[$i]}[0..3]=(split /\t/,$_)[1..4]);
            /snRNA/ && ($snRNA = 1);
            /miRNA/ && ($miRNA = 1);
            $i++;
        }
        my @temp1 = ($SpeType =~ /^[BPV]$/) ? (2..4) : (2..5);
        my @temp2 = ($SpeType =~ /^[BPV]$/) ? (5..7) : (6..9);
        foreach (@temp1) { ${$sample_ncRNA[$_]}[0] .= "(denovo)" }
        foreach (@temp2) { ${$sample_ncRNA[$_]}[0] .= "(homology)" }
#my @index = (-d "$run_dir/$sample/03.ncRNA-Finding/rRNA/homology") ? (1..$#sample_ncRNA) : (($SpeType =~ /^[BPV]$/) ? (1..4,8) : (1..5,10..12) );  #by lss at 20150718
        my @index = (-d "$run_dir/$sample/03.ncRNA-Finding/rRNA/homology") ? (1..$#sample_ncRNA) : (($SpeType =~ /^[BPV]$/) ? (1..4,8) : (1..5,10) );
        $snRNA && (push @index, 11);
        $miRNA && (push @index, 12);
        foreach (@index) { 
            $ncRNA_all .= "$sample\t";
            $ncRNA_all .= join ("\t", @{$sample_ncRNA[$_]});
            $ncRNA_all .= "\n";
        }
    }
	if (-e "$run_dir/$sample/05.Transposon/nuc/$sample.tpsi.nuc.stat"){
		open IN,"$run_dir/$sample/05.Transposon/nuc/$sample.tpsi.nuc.stat"||die$!;
		<IN>;
		while(<IN>){$transposon_all_nuc .= "$_"};
		close IN;
	}	
	if (-e "$run_dir/$sample/05.Transposon/pep/$sample.tpsi.pep.stat"){
		open IN,"$run_dir/$sample/05.Transposon/pep/$sample.tpsi.pep.stat"||die$!;
		<IN>;
		while(<IN>){$transposon_all_pep .= "$_"};
		close IN;
	}	
}
($gene_all) && ($gene_all = "Sample ID\tGenome size(bp)\tGene number(#)\tGene total length(bp)\tGene average length(bp)\tGene length/Genome(%)\n" . $gene_all) && writefile("$outdir/all_sample.gene.stat.xls", $gene_all);
($repbase_all) && ($repbase_all = "Sample ID\tType\tNumber(#)\tTotal Length(bp)\tIn Genome(%)\tAverage length(bp)\n" . $repbase_all) && writefile("$outdir/all_sample.repbase.stat.xls", $repbase_all);
($trf_all) && ($trf_all = "Sample ID\tType\tNumber(#)\tRepeat Size(bp)\tTotal Length(bp)\tIn Genome(%)\n" . $trf_all) && writefile("$outdir/all_sample.trf.stat.xls", $trf_all);
($ncRNA_all) &&  ($ncRNA_all = "Sample ID\tType\tNumber(#)\tAverage length(bp)\tTotal length(bp)\n" . $ncRNA_all) && writefile("$outdir/all_sample.ncRNA.stat.xls", $ncRNA_all);
($transposon_all_nuc) && ( writefile("$outdir/all_sample.transposon.nuc.stat.xls","Sample name\tMatch Number\tMatch length\tAverage length\n$transposon_all_nuc"));
($transposon_all_pep) && ( writefile("$outdir/all_sample.transposon.pep.stat.xls","Sample name\tMatch Number\tMatch length\tAverage length\n$transposon_all_pep"));
#=====sub=============
sub writefile {
    my ($sh_name, $sh) = @_;
    open SH, ">$sh_name" || die "Error: Cannot create $sh_name.\n";
    print SH $sh;
    close SH;
}

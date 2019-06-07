#!/usr/bin/perl

use strict;
use warnings;

if (@ARGV < 2) {
	warn "
	usage: perl delete_dupLine.pl <in seq> <out fseq>\n\n";
	exit;
}

my ($in, $out) = @ARGV;

my %hash_seq;

my $i = 0;
open IN, $in || die;
$/ = "\>";
<IN>;
$/ = "\n";
while (<IN>) {
	chomp;
	next if (/^\s*$/);
	s/^\s+|\s+$//g;
	my $id = $_;
	$/ = "\>";
	chomp (my $seq = <IN>);
	$/ = "\n";
	$hash_seq{$id}{num} = ++$i;
	$hash_seq{$id}{seq} = $seq;
}
close IN;

open OUT, ">$out" || die;
foreach ( sort {$hash_seq{$a}->{num} <=> $hash_seq{$b}->{num}} keys %hash_seq) {
	print OUT ">$_\n$hash_seq{$_}{seq}";
}
close OUT;

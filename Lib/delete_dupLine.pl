#!/usr/bin/perl

use strict;
use warnings;

if (@ARGV < 2) {
	warn "
	usage: perl delete_dupLine.pl <in seq> <out fseq> [pattern]
		line which includes pattern will be saved though it is duplication\n\n";
	exit;
}

my ($in, $out, $pattern) = @ARGV;

my %exist_lines;

open OUT, ">$out" || die;
open IN, $in || die;
while (<IN>) {
	chomp;
	next if (/^\s*$/);
	s/^\s+|\s+$//g;
	if ($pattern && /$pattern/) {
		print OUT "$_\n";
		next;
	}
	if (++$exist_lines{$_} <= 1) {
		print OUT "$_\n";
	}
}
close IN;
close OUT;

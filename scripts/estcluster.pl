#! /usr/bin/perl

use warnings;
use strict;

open (OPEN, "est_cluster.txt");

open (NAME, ">value.txt");
my @array;
while (<OPEN>) {
	chomp $_;
	print $_;
	push (@array, $_);
	}


print $array[5];

close OPEN;
close NAME;

exit;
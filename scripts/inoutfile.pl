#! /usr/bin/perl

use warnings;
use strict;

open (NAME, "$ARGV[0]");

open (OUT, ">nameout.txt");

while (<NAME>) {chomp $_; print OUT "Hello $_!\n";}



close NAME;
close OUT;

exit;
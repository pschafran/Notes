#! /usr/bin/perl

use warnings;
use strict;

open(YO, ">>YO.txt");

my $name = $ARGV[0];

print YO $ARGV[0];

exit;
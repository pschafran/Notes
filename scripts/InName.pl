#! /usr/bin/perl

use warnings;
use strict;

open (IN,"name.txt");
my $name = <IN>;
close IN;

print "Hello $name!\n"; 
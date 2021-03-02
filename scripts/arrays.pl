#! /usr/bin/perl

use warnings;
use strict;

my @genenames;

@genenames = ("ovo", "heartless", "clown", "vasa");

print "The gene names are: @genenames \n";

print "The first gene is: $genenames[0] \n";

my $last_index=$#genenames;

my $num_of_elements = @genenames;

print "Last index: $last_index, Number of elements: $num_of_elements\n";

exit;
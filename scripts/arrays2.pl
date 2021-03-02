#! /usr/bin/perl

use warnings;
use strict;

my @genenames;

@genenames = ("ovo","heartless","clown","vasa");

push (@genenames, "nano");

unshift (@genenames, "runs");

print "The gene list is @genenames\n";

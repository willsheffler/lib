#!/usr/bin/perl -w

use strict;

# Goofy little program that tries to construct a makefile based on scons output.
# Intended to speed up compilation cycles by using make rather than scons.

my $target           = 'default';
my $makefile_name    = 'Makefile';
my @regexes_for_make = qw/ ^g\+\+ ^ld icpc /;

my @makefile_lines;
while ( my $line = <STDIN> ) {
	chomp $line;

	print $line, "\n";
	if ( grep { $line =~ /$_/ } @regexes_for_make ) {
		push @makefile_lines, $line;
	}
}

open  FILE, ">$makefile_name" or die $!;
print FILE "$target:\n";
foreach my $line (@makefile_lines) {
	print FILE "\t", $line, "\n";	
}

close FILE or die $!;

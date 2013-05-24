#!/usr/bin/env perl

# converts and renumbers PDB files
#
# http://mmtsb.scripps.edu/doc/convpdb.pl.html
# 2000, Michael Feig, Brooks group, TSRI

sub usage {
  printf STDERR "usage:   convpdb.pl [options] [PDBfile]\n";
  printf STDERR "options: [-center] [-translate dx dy dz]\n";
  printf STDERR "         [-sel list] [-exclude list]\n";
  printf STDERR "         [-chain id] [-model num] [-nohetero]\n";
  printf STDERR "         [-selseq abbrev]\n";
  printf STDERR "         [-merge pdbfile]\n";
  printf STDERR "         [-renumber start] [-add value]\n";
  printf STDERR "         [-renumberAcrossChains]\n";
  printf STDERR "         [-match pdbfile]\n";
  printf STDERR "         [-setchain id]\n";
  printf STDERR "         [-readseg] [-chainfromseg]\n";
  printf STDERR "         [-charmm19] [-amber]\n";
  printf STDERR "         [-out charmm19 | charmm22 | amber | generic | generic_noh]\n";
  printf STDERR "         [-segnames]\n";
  printf STDERR "         [-fixcoo]\n";
  printf STDERR "         [-ssbond res1:res2[=res1:res2]]\n";
  printf STDERR "         [-findssbonds] [-nossbond]\n";
  printf STDERR "         [-solvate] [-cutoff value]\n";
  printf STDERR "         [-octahedron] [-cubic]\n";
  exit 1;
}

use vars qw ( $perllibdir );

BEGIN {
#  $perllibdir="$ENV{MMTSBDIR}/perl" if (defined $ENV{MMTSBDIR});
   $perllibdir="/work/balej/perl";	
  ($perllibdir=$0)=~s/[^\/]+$// if (!defined $perllibdir);
}

use lib $perllibdir;
use strict;
use Molecule;

my $renumber;
my $renumberAcrossChains;
my $add;
my $fname="-";
my $inmode="";
my $outmode="GENERIC";
my $center=0;
my $sellist;
my $chain;
my $segnames;
my $matchpdb;
my $mergepdb;
my $setchain;
my $ignoreseg=1;
my $chainfromseg=0;
my $fixcoo=0;
my $dx=0.0;
my $dy=0.0;
my $dz=0.0;
my $selmodel=undef;
my $selseq;
my $solvate;
my $shape;
my $cutoff=9.0;
my $excllist;
my $ignorehet=0;
my $ssbonds=();
my $findssbonds=0;
my $nossbond=0;

while ($#ARGV>=0) {
  if ($ARGV[0] eq "-help" || $ARGV[0] eq "-h") {
    &usage();
  } elsif ($ARGV[0] eq "-renumber") {
    shift @ARGV;
    $renumber=shift @ARGV;
  } elsif ($ARGV[0] eq "-renumberAcrossChains") {
    shift @ARGV;
    $renumberAcrossChains=1;
  } elsif ($ARGV[0] eq "-add") {
    shift @ARGV;
    $add=shift @ARGV;
  } elsif ($ARGV[0] eq "-setchain") {
    shift @ARGV;
    $setchain=shift @ARGV;
  } elsif ($ARGV[0] eq "-readseg") {
    shift @ARGV;
    $ignoreseg=0;
  } elsif ($ARGV[0] eq "-chainfromseg") {
    shift @ARGV;
    $chainfromseg=1;
    $ignoreseg=0;
  } elsif ($ARGV[0] eq "-ssbond") {
    shift @ARGV;
    foreach my $s (split(/=/,shift @ARGV)) {
      my @l=split(/:/,$s);
      my $trec={};
      my ($tc1,$tr1)=($l[0]=~/([A-Za-z]*)([\-0-9]*)/);
      $trec->{chain1}=$tc1;
      $trec->{resnum1}=$tr1;
      my ($tc2,$tr2)=($l[1]=~/([A-Za-z]*)([\-0-9]*)/);
      $trec->{chain2}=(defined $tc2 && $tc2 ne "")?$tc2:$tc1;
      $trec->{resnum2}=$tr2;
      push (@{$ssbonds},$trec);
    }
  } elsif ($ARGV[0] eq "-nossbond") {
    shift @ARGV;
    $nossbond=1;
  } elsif ($ARGV[0] eq "-findssbonds") {
    shift @ARGV;
    $findssbonds=1;
  } elsif ($ARGV[0] eq "-charmm19") {
    $inmode="CHARMM19";
    shift @ARGV;
  } elsif ($ARGV[0] eq "-amber") {
    $inmode="AMBER";
    shift @ARGV;
  } elsif ($ARGV[0] eq "-fixcoo") {
    shift @ARGV;
    $fixcoo=1;
  } elsif ($ARGV[0] eq "-out") {
    shift @ARGV;
    $outmode=uc shift @ARGV;
  } elsif ($ARGV[0] eq "-center") {
    shift @ARGV;
    $center=1;
  } elsif ($ARGV[0] eq "-translate") {
    shift @ARGV;
    $dx=shift @ARGV;
    $dy=shift @ARGV;
    $dz=shift @ARGV;
  } elsif ($ARGV[0] eq "-segnames") {
    shift @ARGV;
    $segnames=1;
  } elsif ($ARGV[0] eq "-sel") {
    shift @ARGV;
    $sellist=&GenUtil::fragListFromOption(shift @ARGV);
  } elsif ($ARGV[0] eq "-selseq") {
    shift @ARGV;
    $selseq=shift @ARGV;
  } elsif ($ARGV[0] eq "-exclude") {
    shift @ARGV;
    $excllist=&GenUtil::fragListFromOption(shift @ARGV);
  } elsif ($ARGV[0] eq "-nohetero") {
    shift @ARGV;
    $ignorehet=1;
  } elsif ($ARGV[0] eq "-chain") {
    shift @ARGV;
    $chain=shift @ARGV;
  } elsif ($ARGV[0] eq "-model") {
    shift @ARGV;
    $selmodel=shift @ARGV;
  } elsif ($ARGV[0] eq "-match") {
    shift @ARGV;
    $matchpdb=shift @ARGV;
  } elsif ($ARGV[0] eq "-merge") {
    shift @ARGV;
    $mergepdb=shift @ARGV;
  } elsif ($ARGV[0] eq "-solvate") {
    shift @ARGV;
    $solvate=1;
  } elsif ($ARGV[0] eq "-cubic") {
    shift @ARGV;
    $shape="cubic";
  } elsif ($ARGV[0] eq "-octahedron") {
    shift @ARGV;
    $shape="octahedron";
  } elsif ($ARGV[0] eq "-cutoff" ) {
    shift @ARGV;
    $cutoff=shift @ARGV;
  } else {
    $fname = shift @ARGV;
  }
}

my $mol=Molecule->new();
$mol->readPDB($fname,translate=>$inmode,ignoreseg=>$ignoreseg,
	      ignorehet=>$ignorehet,chainfromseg=>$chainfromseg,
              model=>$selmodel);

$mol->setSSBonds($ssbonds);
$mol->findSSBonds() if ($findssbonds);

my $tseq=Sequence::new($selseq);

$mol->selectChain($chain) if (defined $chain);

$mol->setValidResidues($sellist) if (defined $sellist);
$mol->setValidResidues($excllist,1,1) if (defined $excllist); 

if (defined $selseq) {
  my $foundany=$mol->setValidSequence($selseq);
  die "sequence not found" unless ($foundany);
}
$mol=$mol->clone(1) if (defined $sellist || defined $excllist || defined $selseq || defined $chain);

$mol->setChain($setchain) if (defined $setchain);

if (defined $matchpdb) {
  my $refmol=Molecule->new($matchpdb);
  $mol->match($refmol);
}

if (defined $mergepdb) {
  my $refmol=Molecule->new($mergepdb);
  $mol->merge($refmol);
}

$mol->center() if ($center);
$mol->move($dx,$dy,$dz);
$mol->renumber($renumber) if (defined $renumber);
$mol->renumberAcrossChains() if (defined $renumberAcrossChains);
$mol->shiftResNumber($add)  if (defined $add);
$mol->generateSegNames() if (defined $segnames);
$mol->fixCOO() if ($fixcoo);
$mol->solvate($cutoff,$shape) if (defined $solvate && $solvate);
$mol->writePDB("-",translate=>$outmode,ssbond=>!$nossbond);

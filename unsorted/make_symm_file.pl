#!/usr/bin/perl
##
##
###############################################################################

use strict;
use Math::Trig;   # inv trig ops
use POSIX qw(ceil floor fmod fabs);
#use Getopt::Long qw(permute);
use constant PI    => 4 * atan2(1, 1);

###############################################################################

if ($#ARGV < 0) {
	print STDERR "usage: $0 [options]\n";
	print STDERR "example:   $0 -r 12.0 -c 42.4 41.2 88.6 90.0 90.0 90.0 -s P 1 21 1 -p mystructure.pdb\n";
	print STDERR "options: \n";
	print STDERR "    -r <real>   : [default 8.0] the max CA-CA distance between two interacting chains\n";
	print STDERR "    -b <real>   : Approximate the protein as a ball of this radius (only if no '-p'!)\n";
	print STDERR "    -c <real>x6 : override the unit cell parameters in the PDB with these values\n";
	print STDERR "    -s <string> : override the spacegroup in the PDB with these values\n";
	print STDERR "    -p <string> : Input PDB file (one of -b or -p _must_ be given)\n";
	exit -1;
}

my $pdbfile;
my $interact_dist = 8.0;  # min interaction distance
my $sphere_size = -1; # approx protein as a sphere this size
my @cell_new;
my $spacegp_new;

## parse options (do this by hand since Getopt does not handle this well)
my @suboptions = split /(-[a-z|A-Z] )/, (join ' ',@ARGV);
for ( my $i=0; $i<$#suboptions; $i++ ) {
	if ($suboptions[$i] eq "-r " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		$interact_dist = int( $suboptions[++$i] );
	} elsif ($suboptions[$i] eq "-b " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		$sphere_size = $suboptions[++$i];
	} elsif ($suboptions[$i] eq "-c " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		@cell_new = split /[, ]/,$suboptions[++$i];
	} elsif ($suboptions[$i] eq "-s " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		$spacegp_new = $suboptions[++$i];
	} elsif ($suboptions[$i] eq "-p " && defined $suboptions[$i+1] && $suboptions[$i+1] !~ /^-[a-z|A-Z]/) {
		$pdbfile = $suboptions[++$i];
	} elsif (length($suboptions[$i]) > 0) {
		print STDERR "Error parsing command line while at '".$suboptions[$i]."'\n";
	}
}

#GetOptions('r=f' => \$interact_dist, 'c=f@' => \@cell_new, 's=s@' => \@spacegp_new, 'p=s' => \$pdbfile);

###############################################################################

# crystinfo
my $spacegp = "P 1";
my ($gpid, $nsymm, $Rs, $Ts, $Cs, $cheshire);
my ($A, $B, $C, $alpha, $beta, $gamma) = (0,0,0,90,90,90);
my ($f2c,$c2f);

my @chaintrace;
my @filebuf;
my $CoM = [0,0,0];

if ( length ($pdbfile) > 0 ) {
	###
	### Input PDB file
	open (PDB, $pdbfile) || die "Cannot open $pdbfile.";
	while (<PDB>) {
		chomp;
		if (/^CRYST1/) {
			$A = substr ($_,  6, 9);
			$B = substr ($_, 15, 9);
			$C = substr ($_, 24, 9);
			$alpha = substr ($_, 33, 7);
			$beta  = substr ($_, 40, 7);
			$gamma = substr ($_, 47, 7);
			$spacegp =  substr ($_, 55, 11);
		} elsif (/^ATOM/ || /^HETATM/) {
			my $atom = substr ($_, 12, 4);
			if ($atom eq " CA ") {
				push @chaintrace, [substr ($_, 30, 8),substr ($_, 38, 8),substr ($_, 46, 8)];
				$CoM = vadd( $chaintrace[ $#chaintrace ], $CoM );
			}
	
			push @filebuf, $_;
		}
	}
	close (PDB);

	# override crystal params from cmd line
	if ($#cell_new >= 5) {
		($A,$B,$C,$alpha,$beta,$gamma) = @cell_new;
	}
	if (length($spacegp_new) > 0) {
		$spacegp = $spacegp_new;
	}
	
	# initialize symmops
	($gpid, $nsymm, $Rs, $Ts, $Cs, $cheshire) = spacegp_lookup( $spacegp );
	($f2c,$c2f) = crystparams($A,$B,$C,$alpha,$beta,$gamma);
} else {
	###
	### Interaction radius only
	if ($sphere_size < 0) {
		print STDERR "Must provide an input pdb file with -p _or_ a nonnegative approximate radius with -b!\n";
		exit -1;
	} else {
		print STDERR "Warning! No input structure provided ... generating symmetry mates using interaction radius only!\n";
	}

	# override crystal params from cmd line
	if ($#cell_new >= 5) {
		($A,$B,$C,$alpha,$beta,$gamma) = @cell_new;
	}
	if (length($spacegp_new) > 0) {
		$spacegp = $spacegp_new;
	}
	
	# initialize symmops
	($gpid, $nsymm, $Rs, $Ts, $cheshire) = spacegp_lookup( $spacegp );
	($f2c,$c2f) = crystparams($A,$B,$C,$alpha,$beta,$gamma);

	# fake atom
	my $fake_fX = [ ($cheshire->[0][0]+$cheshire->[0][1])/2 ,
	                ($cheshire->[1][0]+$cheshire->[1][1])/2 ,
	                ($cheshire->[2][0]+$cheshire->[2][1])/2 ];
	my $fake_cX = mapply( $f2c , $fake_fX );
	push @chaintrace, $fake_cX;
	my $fakePDBline = sprintf "ATOM      1  CA  ALA A   1     %7.3f %7.3f %7.3f  1.00  0.00           C  ", 
	                          $fake_cX->[0], $fake_cX->[1], $fake_cX->[2];
	push @filebuf, $fakePDBline;
	$CoM = $fake_cX;

	if ($sphere_size > 0) {
		my $plusZed = [0,0,$sphere_size];
		my $nBetaSteps = int( PI*$sphere_size / 8.0 )+1;
		my $counter = 2;
		foreach my $i (1..$nBetaSteps) {
			my $beta = PI*($i-0.5)/($nBetaSteps);
			my $nGammaSteps = int( PI*$sphere_size*sin( $beta )/ 4.0 )+1;

			foreach my $j (1..$nGammaSteps) {
				my $gamma = 2*PI*($j-0.5)/($nGammaSteps);

				my $rot_ij = euler( 0,$beta,$gamma );
				my $cX_ij = vadd( mapply( $rot_ij , $plusZed ) , $fake_cX );
				push @chaintrace, $cX_ij;
				$fakePDBline = sprintf "ATOM    %3d  H   ALA A   1     %7.3f %7.3f %7.3f  1.00  0.00           H  ", 
										  $counter++, $cX_ij->[0], $cX_ij->[1], $cX_ij->[2];
				push @filebuf, $fakePDBline;
			}
		}
	}
}

my $nAtms = $#chaintrace+1;
$CoM = [ $CoM->[0]/$nAtms , $CoM->[1]/$nAtms , $CoM->[2]/$nAtms ];

# find residue closest to CoM of the system
my $minDist2 = 9999999;
my $minRes = 1;
foreach my $i (0..$nAtms) {
	my $dist2 = vnorm2( vsub( $CoM, $chaintrace[ $i ] ) );
	if ($dist2 < $minDist2) {
		$minDist2 = $dist2;
		$minRes = $i+1;
	}
}

# frac coords of chaintrace
my @fchaintrace;
foreach my $X_i (@chaintrace) {
	my $fX_i = mapply( $c2f , $X_i );
	push @fchaintrace , $fX_i;
}


# do the symmetric expansion
my %symminterface = ();
foreach my $fX_i (@fchaintrace) {
	foreach my $fY_i (@fchaintrace) {
		foreach my $j_symm (0..($nsymm-1)) {
			my $fY_j = vadd( mapply($Rs->[$j_symm],$fY_i) , $Ts->[$j_symm] );
			my $delfXY = vsub( $fY_j,$fX_i );
			my $delfXY_min = vminmod( $delfXY , [1.0,1.0,1.0] );

			# distance check ...
			my $delXY = mapply( $f2c, $delfXY_min );
			my $dist2XY = vnorm2( $delXY );

			if ($dist2XY < $interact_dist*$interact_dist) {
				my $shiftXY = vsub( $delfXY_min , $delfXY );  # should be all integers

				# we have a hit! tag symmop $jsymm and offset $shiftXY as a symmetic interface!
				my $id = $j_symm.",".floor($shiftXY->[0]+0.5).",".floor($shiftXY->[1]+0.5).",".floor($shiftXY->[2]+0.5);
				if (!defined $symminterface{ $id }) {
					$symminterface{ $id } = $j_symm+vnorm2($shiftXY);
				}
			}
		}
	}
}

# symmetric interfaces will always come in pairs
# find these pairs and place in syminterfaces_paired like [0, 1, 1', 2, 2', ...]
my @syminterfaces = sort { $symminterface{$a} <=> $symminterface{$b} } keys %symminterface;
my @syminterfaces_paired;
my @syminterfaces_unpaired;
push @syminterfaces_unpaired, $syminterfaces[0];   # first the origin
my %paired = (0=>1);

OUTER: foreach my $i (1..$#syminterfaces-1) {
	next if (defined $paired{$i});

	my ($symm_i,$shiftX_i,$shiftY_i,$shiftZ_i) = split ',',$syminterfaces[ $i ];
	my $R_i = $Rs->[$symm_i];
	my $T_i =vadd( $Ts->[$symm_i], [$shiftX_i,$shiftY_i,$shiftZ_i] );

	foreach my $j ($i+1..$#syminterfaces) {
		my ($symm_j,$shiftX_j,$shiftY_j,$shiftZ_j) = split ',',$syminterfaces[ $j ];
		my $R_j = $Rs->[$symm_j];
		my $T_j =vadd( $Ts->[$symm_j], [$shiftX_j,$shiftY_j,$shiftZ_j] );

		# if transform i is the inverse of transform j we have our pair
		if ( is_inverse( $R_i,$T_i, $R_j,$T_j ) ) {
			$paired{ $j } = 1;
			push @syminterfaces_paired, $syminterfaces[$i];
			push @syminterfaces_paired, $syminterfaces[$j];
			next OUTER;
		}
	}
	#print STDERR "Interface $i (".$syminterfaces[ $i ].") unpaired!\n";
	push @syminterfaces_unpaired, $syminterfaces[$i];   # self-paired
}

#######################################
##
## write output symm file
# symmetry_name c4
# subunits 4
# number_of_interfaces 2
# E = 3*E2
# anchor_residue 17
my $symmname = $pdbfile;
$symmname =~ s/\.pdb$//;
$symmname = $symmname." $spacegp";
$symmname =~ s/\s*$//;
$symmname =~ s/\s/_/g;
print "symmetry_name $symmname\n";

my $nsymm = $#syminterfaces_paired+1+$#syminterfaces_unpaired+1;
print "subunits $nsymm\n";

my $ninterfaces = (1+$#syminterfaces_paired)/2 + $#syminterfaces_unpaired;
print "number_of_interfaces $ninterfaces\n";

my $firstterm = 1;
print "E = ";
foreach my $i (2..(1+$#syminterfaces_unpaired)) {
	if ($firstterm == 0) {
		print " + ";
	} else {
		$firstterm = 0;
	}
	print "E".($i);
}
foreach my $i (1..(1+$#syminterfaces_paired)/2) {
	if ($firstterm == 0) {
		print " + ";
	} else {
		$firstterm = 0;
	}
	print "2*E".(2*$i+$#syminterfaces_paired);
}
print "\n";
print "anchor_residue $minRes\n";

# virtual_coordinates_start
# xyz -1,0,0 0,1,0 0,0,0
# xyz 0,-1,0 -1,0,0 0,0,0
# xyz 1,0,0 0,-1,0 0,0,0
# xyz 0,1,0 1,0,0 0,0,0
# virtual_coordinates_stop
print "virtual_coordinates_start\n";
my @syminterfaces_all;
push @syminterfaces_all, @syminterfaces_unpaired;
push @syminterfaces_all, @syminterfaces_paired;
foreach my $symmkey (@syminterfaces_all) {
	my ($j_symm,$shiftX,$shiftY,$shiftZ) = split ',',$symmkey;
	print "xyz";

	# X
	my $fX  = mapply( $c2f , [1,0,0] );
	my $sfX  = mapply($Rs->[$j_symm],$fX);
	my $sX = mapply( $f2c , $sfX );
	my $string = sprintf("%.3f,%.3f,%.3f", $sX->[0], $sX->[1], $sX->[2]);
	print " ".$string;

	# Y
	$fX  = mapply( $c2f , [0,1,0] );
	$sfX  = mapply($Rs->[$j_symm],$fX);
	$sX = mapply( $f2c , $sfX );
	my $string = sprintf("%.3f,%.3f,%.3f", $sX->[0], $sX->[1], $sX->[2]);
	print " ".$string;

	# orig
	$sfX  = vadd( mapply($Rs->[$j_symm],[0,0,0]) , $Ts->[$j_symm] );
	$sfX = vadd( $sfX , [$shiftX,$shiftY,$shiftZ] );
	$sX = mapply( $f2c , $sfX );

	my $string = sprintf("%.3f,%.3f,%.3f", $sX->[0], $sX->[1], $sX->[2]);
	print " ".$string;
	print "\n";
}
print "virtual_coordinates_stop\n";

# connect_virtual 1 2
# connect_virtual 3 4
# connect_virtual 1 3
foreach my $i (2..$nsymm) {
	print "connect_virtual 1 $i\n";
}

# set_dof 1 x angle_x angle_y angle_z
print "set_dof 1";
if ($cheshire->[0][0] != $cheshire->[0][1]) { print " x"; }
if ($cheshire->[1][0] != $cheshire->[1][1]) { print " y"; }
if ($cheshire->[2][0] != $cheshire->[2][1]) { print " z"; }
print " angle_x angle_y angle_z\n";


########################################
##
## write output pdb
my $outpdb = $pdbfile;
if ($outpdb =~ /\.pdb$/) {
	$outpdb =~ s/\.pdb$/_symm.pdb/;
} else {
	$outpdb = $outpdb."_symm.pdb";
}
open (OUTPDB, ">$outpdb");

my $chnidx = 0;
my $chains = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
foreach my $symmkey (@syminterfaces_all) {
	my ($j_symm,$shiftX,$shiftY,$shiftZ) = split ',',$symmkey;
	foreach my $line (@filebuf) {
		my $linecopy = $line;

		my $X = [substr ($line, 30, 8),substr ($line, 38, 8),substr ($line, 46, 8)];
		my $fX = mapply( $c2f , $X );
		my $sfX  = vadd( mapply($Rs->[$j_symm],$fX) , $Ts->[$j_symm] );
		my $sfX1 = vadd( $sfX , [$shiftX,$shiftY,$shiftZ] );
		my $sX = mapply( $f2c , $sfX1 );

		substr ($linecopy, 30, 8) = sprintf ("%8.3f", $sX->[0]);
		substr ($linecopy, 38, 8) = sprintf ("%8.3f", $sX->[1]);
		substr ($linecopy, 46, 8) = sprintf ("%8.3f", $sX->[2]);
		substr ($linecopy, 21, 1) = substr ($chains, $chnidx, 1);

		print OUTPDB $linecopy."\n";
	}
	print OUTPDB "TER   \n";
	$chnidx++;
}
close(OUTPDB);


exit 0;

########
# subs #
########

# rotation from euler angles
sub euler {
	my ($aa, $bb, $gg) = @_;
	my $MM;

	$MM->[0][0] = (-sin($aa)*cos($bb)*sin($gg) + cos($aa)*cos($gg));
	$MM->[0][1] = ( cos($aa)*cos($bb)*sin($gg) + sin($aa)*cos($gg));
	$MM->[0][2] = ( sin($bb)*sin($gg));
	$MM->[1][0] = (-sin($aa)*cos($bb)*cos($gg) - cos($aa)*sin($gg));
	$MM->[1][1] = ( cos($aa)*cos($bb)*cos($gg) - sin($aa)*sin($gg));
	$MM->[1][2] = ( sin($bb)*cos($gg));
	$MM->[2][0] = ( sin($aa)*sin($bb));
	$MM->[2][1] = (-cos($aa)*sin($bb));
	$MM->[2][2] = ( cos($bb));

	return $MM;
}

# vector addition
sub vadd {
	my ($x, $y) = @_;
	return [ $x->[0]+$y->[0], $x->[1]+$y->[1], $x->[2]+$y->[2] ];
}

# vector subtraction
sub vsub {
	my ($x, $y) = @_;
	return [ $x->[0]-$y->[0], $x->[1]-$y->[1], $x->[2]-$y->[2] ];
}

# "min mod"
sub minmod {
	my ($x,$y) = @_;
	my $r = fmod($x,$y);
	if ($r < -fabs( $y/2.0 ) ) { $r += fabs( $y ); }
	elsif ($r >  fabs( $y/2.0 ) ) { $r -= fabs( $y ); }
	return $r;
}

# vector min-modulus
sub vminmod {
	my ($x,$y) = @_;
	return [ minmod($x->[0],$y->[0]), minmod($x->[1],$y->[1]), minmod($x->[2],$y->[2]) ];
}

# matrix x vector mult
sub mapply {
	my ($rotmat, $cart) = @_;
	my $out = [0, 0, 0];
	my ($i, $j);
	for ($i=0; $i < 3; ++$i) {
		for ($j=0; $j < 3; ++$j) {
			$out->[$i] += $rotmat->[$i][$j] * $cart->[$j];
		}
	}
	return $out;
}

# matrix x matrix mult
sub mmult {
	my ($m1, $m2) = @_;
	my $out = [ [0,0,0], [0,0,0], [0,0,0] ];
	my ($i, $j, $k);
	for ($i=0; $i<3; ++$i) {
		for ($j=0; $j<3; ++$j) {
			for ($k=0; $k<3; ++$k) {
				$out->[$i][$j] += $m1->[$i][$k] * $m2->[$k][$j];
			}
		}
	}
	return $out;
}



# matrix inversion
sub minv {
	my $M = shift;
	my $Minv = [ [1,0,0] , [0,1,0] , [0,0,1] ];
	my $D = $M->[0][0] * ( $M->[1][1]*$M->[2][2] - $M->[2][1]*$M->[1][2] ) -
		    $M->[0][1] * ( $M->[1][0]*$M->[2][2] - $M->[1][2]*$M->[2][0] ) +
		    $M->[0][2] * ( $M->[1][0]*$M->[2][1] - $M->[1][1]*$M->[2][0] );
	if ($D == 0)  {
		print STDERR "ERROR ... Inversion of singular matrix!\n";
		exit -1;
	}

	$Minv->[0][0] =  ($M->[1][1]*$M->[2][2]-$M->[1][2]*$M->[2][1])/$D;
	$Minv->[0][1] = -($M->[0][1]*$M->[2][2]-$M->[0][2]*$M->[2][1])/$D;
	$Minv->[0][2] =  ($M->[0][1]*$M->[1][2]-$M->[0][2]*$M->[1][1])/$D;
	$Minv->[1][0] = -($M->[1][0]*$M->[2][2]-$M->[2][0]*$M->[1][2])/$D;
	$Minv->[1][1] =  ($M->[0][0]*$M->[2][2]-$M->[0][2]*$M->[2][0])/$D;
	$Minv->[1][2] = -($M->[0][0]*$M->[1][2]-$M->[0][2]*$M->[1][0])/$D;
	$Minv->[2][0] =  ($M->[1][0]*$M->[2][1]-$M->[2][0]*$M->[1][1])/$D;
	$Minv->[2][1] = -($M->[0][0]*$M->[2][1]-$M->[0][1]*$M->[2][0])/$D;
	$Minv->[2][2] =  ($M->[0][0]*$M->[1][1]-$M->[0][1]*$M->[1][0])/$D;

	return $Minv;
}

# vector norm
sub vnorm {
	my $x = shift;
	return sqrt( ($x->[0]*$x->[0]) + ($x->[1]*$x->[1]) + ($x->[2]*$x->[2]) );
}

# vector norm^2
sub vnorm2 {
	my $x = shift;
	return ( ($x->[0]*$x->[0]) + ($x->[1]*$x->[1]) + ($x->[2]*$x->[2]) );
}

# cart distance
sub vdist {
	my ($x1, $x2) = @_;
	return sqrt (
	           ($x1->[0]-$x2->[0])*($x1->[0]-$x2->[0]) +
	           ($x1->[1]-$x2->[1])*($x1->[1]-$x2->[1]) +
	           ($x1->[2]-$x2->[2])*($x1->[2]-$x2->[2])
	            );
}

# are two transformations the inverso of one another
sub is_inverse {
	my $tol = 1e-8;
	my ( $R_i,$T_i, $R_j,$T_j ) = @_;

	my $testR = mmult( $R_i , $R_j );
	my $testT = vadd( mapply( $R_j,$T_i ) , $T_j );

	my $errR = square($testR->[0][0]-1) + square($testR->[1][1]-1) + square($testR->[2][2]-1) +
	           square($testR->[0][1])   + square($testR->[0][2])   + square($testR->[1][2]) +
	           square($testR->[1][0])   + square($testR->[2][0])   + square($testR->[2][1]);
	my $errT = square($testT->[0])+square($testT->[1])+square($testT->[2]);

	if ($errR < $tol && $errT < $tol) { return 1; }
	return 0;
}

###########
# f2c/c2f #
###########

sub d2r { return (@_[0]*PI/180); }
sub square { return @_[0]*@_[0]; }

# my ($f2c,$c2f) = crystparams($a,$b,$c,$alpha,$beta,$gamma)
sub crystparams {
	my ($a,$b,$c,$alpha,$beta,$gamma) = @_;

	if ($a*$b*$c == 0) {
		print STDERR "Must provide valid crystal parameters!\n";
		exit -1;
	}

	my $f2c = [ [0,0,0] , [0,0,0] , [0,0,0] ];

	my $ca = cos(d2r($alpha)); my $cb = cos(d2r($beta)); my $cg = cos(d2r($gamma));
	my $sa = sin(d2r($alpha)); my $sb = sin(d2r($beta)); my $sg = sin(d2r($gamma));
	$f2c->[0][0] = $A;
	$f2c->[0][1] = $B * $cg;
	$f2c->[0][2] = $C * $cb;
	$f2c->[1][1] = $B * $sg;
	$f2c->[1][2] = $C * ($ca - $cb*$cg) / $sg;
	$f2c->[2][2] = $C * $sb * sqrt(1.0 - square(($cb*$cg - $ca)/($sb*$sg)));

	my $c2f = minv($f2c);

	return ($f2c,$c2f);
}

#################
# symmgp lookup #
#################

# my ($gpid, $nsymm, $Rs, $Ts, $cheshire) = spacegp( $spacegpid );
sub spacegp_lookup {
	my $id = shift @_;
	my ($gpid, $nsymm, @Rs, @Ts, @Cs, @R_all, @T_all, $cheshire);

	# strip leading/trailing spaces
	$id =~ s/^\s+//;
	$id =~ s/\s+$//;
	
	if ($id eq "P 1") {
		$gpid = 1;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,0] , [0,0] , [0,0] ];
	}
	elsif ($id eq "P -1") {
		$gpid = 2;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 1 2 1") {
		$gpid = 3;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,0] , [0,1/2] ];
	}
	elsif ($id eq "P 1 1 2") {
		$gpid = 3;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,0] , [0,1/2] ];
	}
	elsif ($id eq "P 2 1 1") {
		$gpid = 3;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,0] , [0,1/2] ];
	}
	elsif ($id eq "P 1 21 1") {
		$gpid = 4;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,0] , [0,1/2] ];
	}
	elsif ($id eq "P 1 1 21") {
		$gpid = 4;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,0] , [0,1/2] ];
	}
	elsif ($id eq "P 21 1 1") {
		$gpid = 4;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,0] , [0,1/2] ];
	}
	elsif ($id eq "C 1 2 1") {
		$gpid = 5;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,0] , [0,1/2] ];
	}
	elsif ($id eq "A 1 2 1") {
		$gpid = 5;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,0] , [0,1/2] ];
	}
	elsif ($id eq "I 1 2 1") {
		$gpid = 5;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,0] , [0,1/2] ];
	}
	elsif ($id eq "A 1 1 2") {
		$gpid = 5;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,0] , [0,1/2] ];
	}
	elsif ($id eq "B 1 1 2") {
		$gpid = 5;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,0] , [0,1/2] ];
	}
	elsif ($id eq "I 1 1 2") {
		$gpid = 5;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,0] , [0,1/2] ];
	}
	elsif ($id eq "B 2 1 1") {
		$gpid = 5;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,0] , [0,1/2] ];
	}
	elsif ($id eq "C 2 1 1") {
		$gpid = 5;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,0] , [0,1/2] ];
	}
	elsif ($id eq "I 2 1 1") {
		$gpid = 5;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,0] , [0,1/2] ];
	}
	elsif ($id eq "P 1 m 1") {
		$gpid = 6;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 1 1 m") {
		$gpid = 6;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P m 1 1") {
		$gpid = 6;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 1 c 1") {
		$gpid = 7;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 1 n 1") {
		$gpid = 7;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 1 a 1") {
		$gpid = 7;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 1 1 a") {
		$gpid = 7;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 1 1 n") {
		$gpid = 7;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 1 1 b") {
		$gpid = 7;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P b 1 1") {
		$gpid = 7;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P n 1 1") {
		$gpid = 7;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P c 1 1") {
		$gpid = 7;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "C 1 m 1") {
		$gpid = 8;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "A 1 m 1") {
		$gpid = 8;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I 1 m 1") {
		$gpid = 8;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "A 1 1 m") {
		$gpid = 8;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "B 1 1 m") {
		$gpid = 8;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I 1 1 m") {
		$gpid = 8;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "B m 1 1") {
		$gpid = 8;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "C m 1 1") {
		$gpid = 8;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I m 1 1") {
		$gpid = 8;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "C 1 c 1") {
		$gpid = 9;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "A 1 n 1") {
		$gpid = 9;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I 1 a 1") {
		$gpid = 9;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "A 1 a 1") {
		$gpid = 9;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "C 1 n 1") {
		$gpid = 9;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I 1 c 1") {
		$gpid = 9;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "A 1 1 a") {
		$gpid = 9;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "B 1 1 n") {
		$gpid = 9;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I 1 1 b") {
		$gpid = 9;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "B 1 1 b") {
		$gpid = 9;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "A 1 1 n") {
		$gpid = 9;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I 1 1 a") {
		$gpid = 9;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "B b 1 1") {
		$gpid = 9;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "C n 1 1") {
		$gpid = 9;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I c 1 1") {
		$gpid = 9;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "C c 1 1") {
		$gpid = 9;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "B n 1 1") {
		$gpid = 9;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I b 1 1") {
		$gpid = 9;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,0] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 1 2/m 1") {
		$gpid = 10;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 1 1 2/m") {
		$gpid = 10;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 2/m 1 1") {
		$gpid = 10;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 1 21/m 1") {
		$gpid = 11;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 1 1 21/m") {
		$gpid = 11;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 21/m 1 1") {
		$gpid = 11;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "C 1 2/m 1") {
		$gpid = 12;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "A 1 2/m 1") {
		$gpid = 12;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 1 2/m 1") {
		$gpid = 12;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "A 1 1 2/m") {
		$gpid = 12;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "B 1 1 2/m") {
		$gpid = 12;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 1 1 2/m") {
		$gpid = 12;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "B 2/m 1 1") {
		$gpid = 12;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "C 2/m 1 1") {
		$gpid = 12;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 2/m 1 1") {
		$gpid = 12;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 1 2/c 1") {
		$gpid = 13;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 1 2/n 1") {
		$gpid = 13;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 1 2/a 1") {
		$gpid = 13;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 1 1 2/a") {
		$gpid = 13;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 1 1 2/n") {
		$gpid = 13;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 1 1 2/b") {
		$gpid = 13;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 2/b 1 1") {
		$gpid = 13;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 2/n 1 1") {
		$gpid = 13;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 2/c 1 1") {
		$gpid = 13;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 1 21/c 1") {
		$gpid = 14;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 1 21/n 1") {
		$gpid = 14;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 1 21/a 1") {
		$gpid = 14;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 1 1 21/a") {
		$gpid = 14;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 1 1 21/n") {
		$gpid = 14;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 1 1 21/b") {
		$gpid = 14;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 21/b 1 1") {
		$gpid = 14;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 21/n 1 1") {
		$gpid = 14;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 21/c 1 1") {
		$gpid = 14;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "C 1 2/c 1") {
		$gpid = 15;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "A 1 2/n 1") {
		$gpid = 15;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 1 2/a 1") {
		$gpid = 15;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "A 1 2/a 1") {
		$gpid = 15;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "C 1 2/n 1") {
		$gpid = 15;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 1 2/c 1") {
		$gpid = 15;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "A 1 1 2/a") {
		$gpid = 15;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "B 1 1 2/n") {
		$gpid = 15;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 1 1 2/b") {
		$gpid = 15;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "B 1 1 2/b") {
		$gpid = 15;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "A 1 1 2/n") {
		$gpid = 15;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 1 1 2/a") {
		$gpid = 15;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "B 2/b 1 1") {
		$gpid = 15;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "C 2/n 1 1") {
		$gpid = 15;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 2/c 1 1") {
		$gpid = 15;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "C 2/c 1 1") {
		$gpid = 15;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "B 2/n 1 1") {
		$gpid = 15;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 2/b 1 1") {
		$gpid = 15;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 2 2 2") {
		$gpid = 16;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 2 2 21") {
		$gpid = 17;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 21 2 2") {
		$gpid = 17;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 2 21 2") {
		$gpid = 17;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 21 21 2") {
		$gpid = 18;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 2 21 21") {
		$gpid = 18;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 21 2 21") {
		$gpid = 18;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 21 21 21") {
		$gpid = 19;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "C 2 2 21") {
		$gpid = 20;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "A 21 2 2") {
		$gpid = 20;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "B 2 21 2") {
		$gpid = 20;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "C 2 2 2") {
		$gpid = 21;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "A 2 2 2") {
		$gpid = 21;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "B 2 2 2") {
		$gpid = 21;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "F 2 2 2") {
		$gpid = 22;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 2 2 2") {
		$gpid = 23;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 21 21 21") {
		$gpid = 24;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P m m 2") {
		$gpid = 25;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 2 m m") {
		$gpid = 25;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P m 2 m") {
		$gpid = 25;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P m c 21") {
		$gpid = 26;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P c m 21") {
		$gpid = 26;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 21 m a") {
		$gpid = 26;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 21 a m") {
		$gpid = 26;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P b 21 m") {
		$gpid = 26;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P m 21 b") {
		$gpid = 26;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P c c 2") {
		$gpid = 27;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 2 a a") {
		$gpid = 27;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P b 2 b") {
		$gpid = 27;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P m a 2") {
		$gpid = 28;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P b m 2") {
		$gpid = 28;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 2 m b") {
		$gpid = 28;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 2 c m") {
		$gpid = 28;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P c 2 m") {
		$gpid = 28;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P m 2 a") {
		$gpid = 28;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P c a 21") {
		$gpid = 29;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P b c 21") {
		$gpid = 29;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 21 a b") {
		$gpid = 29;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 21 c a") {
		$gpid = 29;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P c 21 b") {
		$gpid = 29;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P b 21 a") {
		$gpid = 29;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P n c 2") {
		$gpid = 30;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P c n 2") {
		$gpid = 30;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 2 n a") {
		$gpid = 30;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 2 a n") {
		$gpid = 30;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P b 2 n") {
		$gpid = 30;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P n 2 b") {
		$gpid = 30;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P m n 21") {
		$gpid = 31;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P n m 21") {
		$gpid = 31;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 21 m n") {
		$gpid = 31;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 21 n m") {
		$gpid = 31;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P n 21 m") {
		$gpid = 31;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P m 21 n") {
		$gpid = 31;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P b a 2") {
		$gpid = 32;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 2 c b") {
		$gpid = 32;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P c 2 a") {
		$gpid = 32;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P n a 21") {
		$gpid = 33;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P b n 21") {
		$gpid = 33;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 21 n b") {
		$gpid = 33;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 21 c n") {
		$gpid = 33;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P c 21 n") {
		$gpid = 33;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P n 21 a") {
		$gpid = 33;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P n n 2") {
		$gpid = 34;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 2 n n") {
		$gpid = 34;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P n 2 n") {
		$gpid = 34;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "C m m 2") {
		$gpid = 35;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "A 2 m m") {
		$gpid = 35;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "B m 2 m") {
		$gpid = 35;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "C m c 21") {
		$gpid = 36;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "C c m 21") {
		$gpid = 36;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "A 21 m a") {
		$gpid = 36;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "A 21 a m") {
		$gpid = 36;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "B b 21 m") {
		$gpid = 36;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "B m 21 b") {
		$gpid = 36;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "C c c 2") {
		$gpid = 37;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "A 2 a a") {
		$gpid = 37;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "B b 2 b") {
		$gpid = 37;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "A m m 2") {
		$gpid = 38;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "B m m 2") {
		$gpid = 38;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "B 2 m m") {
		$gpid = 38;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "C 2 m m") {
		$gpid = 38;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "C m 2 m") {
		$gpid = 38;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "A m 2 m") {
		$gpid = 38;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "A b m 2") {
		$gpid = 39;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "B m a 2") {
		$gpid = 39;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "B 2 c m") {
		$gpid = 39;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "C 2 m b") {
		$gpid = 39;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "C m 2 a") {
		$gpid = 39;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "A c 2 m") {
		$gpid = 39;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "A m a 2") {
		$gpid = 40;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "B b m 2") {
		$gpid = 40;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "B 2 m b") {
		$gpid = 40;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "C 2 c m") {
		$gpid = 40;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "C c 2 m") {
		$gpid = 40;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "A m 2 a") {
		$gpid = 40;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "A b a 2") {
		$gpid = 41;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "B b a 2") {
		$gpid = 41;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "B 2 c b") {
		$gpid = 41;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "C 2 c b") {
		$gpid = 41;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "C c 2 a") {
		$gpid = 41;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "A c 2 a") {
		$gpid = 41;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "F m m 2") {
		$gpid = 42;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "F 2 m m") {
		$gpid = 42;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "F m 2 m") {
		$gpid = 42;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "F d d 2") {
		$gpid = 43;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.75,0.75,0.25];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "F 2 d d") {
		$gpid = 43;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.25,0.75,0.75];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "F d 2 d") {
		$gpid = 43;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.75,0.25,0.75];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I m m 2") {
		$gpid = 44;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I 2 m m") {
		$gpid = 44;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I m 2 m") {
		$gpid = 44;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I b a 2") {
		$gpid = 45;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I 2 c b") {
		$gpid = 45;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I c 2 a") {
		$gpid = 45;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I m a 2") {
		$gpid = 46;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I b m 2") {
		$gpid = 46;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I 2 m b") {
		$gpid = 46;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I 2 c m") {
		$gpid = 46;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I c 2 m") {
		$gpid = 46;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I m 2 a") {
		$gpid = 46;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P m m m") {
		$gpid = 47;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P n n n :1") {
		$gpid = 48;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P n n n :2") {
		$gpid = 48;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P c c m") {
		$gpid = 49;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P m a a") {
		$gpid = 49;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P b m b") {
		$gpid = 49;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P b a n :1") {
		$gpid = 50;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P b a n :2") {
		$gpid = 50;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P n c b :1") {
		$gpid = 50;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P n c b :2") {
		$gpid = 50;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P c n a :1") {
		$gpid = 50;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P c n a :2") {
		$gpid = 50;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P m m a") {
		$gpid = 51;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P m m b") {
		$gpid = 51;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P b m m") {
		$gpid = 51;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P c m m") {
		$gpid = 51;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P m c m") {
		$gpid = 51;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P m a m") {
		$gpid = 51;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P n n a") {
		$gpid = 52;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P n n b") {
		$gpid = 52;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P b n n") {
		$gpid = 52;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P c n n") {
		$gpid = 52;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P n c n") {
		$gpid = 52;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P n a n") {
		$gpid = 52;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P m n a") {
		$gpid = 53;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P n m b") {
		$gpid = 53;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P b m n") {
		$gpid = 53;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P c n m") {
		$gpid = 53;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P n c m") {
		$gpid = 53;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P m a n") {
		$gpid = 53;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P c c a") {
		$gpid = 54;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P c c b") {
		$gpid = 54;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P b a a") {
		$gpid = 54;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P c a a") {
		$gpid = 54;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P b c b") {
		$gpid = 54;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P b a b") {
		$gpid = 54;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P b a m") {
		$gpid = 55;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P m c b") {
		$gpid = 55;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P c m a") {
		$gpid = 55;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P c c n") {
		$gpid = 56;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P n a a") {
		$gpid = 56;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P b n b") {
		$gpid = 56;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P b c m") {
		$gpid = 57;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P c a m") {
		$gpid = 57;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P m c a") {
		$gpid = 57;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P m a b") {
		$gpid = 57;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P b m a") {
		$gpid = 57;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P c m b") {
		$gpid = 57;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P n n m") {
		$gpid = 58;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P m n n") {
		$gpid = 58;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P n m n") {
		$gpid = 58;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P m m n :1") {
		$gpid = 59;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P m m n :2") {
		$gpid = 59;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P n m m :1") {
		$gpid = 59;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P n m m :2") {
		$gpid = 59;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P m n m :1") {
		$gpid = 59;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P m n m :2") {
		$gpid = 59;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P b c n") {
		$gpid = 60;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P c a n") {
		$gpid = 60;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P n c a") {
		$gpid = 60;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P n a b") {
		$gpid = 60;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P b n a") {
		$gpid = 60;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P c n b") {
		$gpid = 60;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P b c a") {
		$gpid = 61;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P c a b") {
		$gpid = 61;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P n m a") {
		$gpid = 62;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P m n b") {
		$gpid = 62;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P b n m") {
		$gpid = 62;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P c m n") {
		$gpid = 62;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P m c n") {
		$gpid = 62;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P n a m") {
		$gpid = 62;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "C m c m") {
		$gpid = 63;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "C c m m") {
		$gpid = 63;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "A m m a") {
		$gpid = 63;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "A m a m") {
		$gpid = 63;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "B b m m") {
		$gpid = 63;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "B m m b") {
		$gpid = 63;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "C m c a") {
		$gpid = 64;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "C c m b") {
		$gpid = 64;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "A b m a") {
		$gpid = 64;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "A c a m") {
		$gpid = 64;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "B b c m") {
		$gpid = 64;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "B m a b") {
		$gpid = 64;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "C m m m") {
		$gpid = 65;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "A m m m") {
		$gpid = 65;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "B m m m") {
		$gpid = 65;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "C c c m") {
		$gpid = 66;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "A m a a") {
		$gpid = 66;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "B b m b") {
		$gpid = 66;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "C m m a") {
		$gpid = 67;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "C m m b") {
		$gpid = 67;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "A b m m") {
		$gpid = 67;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "A c m m") {
		$gpid = 67;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "B m c m") {
		$gpid = 67;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "B m a m") {
		$gpid = 67;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "C c c a :1") {
		$gpid = 68;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "C c c a :2") {
		$gpid = 68;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "C c c b :1") {
		$gpid = 68;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "C c c b :2") {
		$gpid = 68;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "A b a a :1") {
		$gpid = 68;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "A b a a :2") {
		$gpid = 68;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "A c a a :1") {
		$gpid = 68;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "A c a a :2") {
		$gpid = 68;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "B b c b :1") {
		$gpid = 68;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "B b c b :2") {
		$gpid = 68;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "B b a b :1") {
		$gpid = 68;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "B b a b :2") {
		$gpid = 68;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "F m m m") {
		$gpid = 69;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "F d d d :1") {
		$gpid = 70;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.75,0.75,0.25];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.75,0.25,0.75];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "F d d d :2") {
		$gpid = 70;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.25,0.25,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.25,0.25];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.25,0,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.75,0.75,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.75,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.75,0,0.75];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I m m m") {
		$gpid = 71;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I b a m") {
		$gpid = 72;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I m c b") {
		$gpid = 72;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I c m a") {
		$gpid = 72;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I b c a") {
		$gpid = 73;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I c a b") {
		$gpid = 73;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I m m a") {
		$gpid = 74;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I m m b") {
		$gpid = 74;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I b m m") {
		$gpid = 74;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I c m m") {
		$gpid = 74;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I m c m") {
		$gpid = 74;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I m a m") {
		$gpid = 74;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 4") {
		$gpid = 75;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 41") {
		$gpid = 76;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.75];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 42") {
		$gpid = 77;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 43") {
		$gpid = 78;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.75];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.25];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I 4") {
		$gpid = 79;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I 41") {
		$gpid = 80;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.75];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P -4") {
		$gpid = 81;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I -4") {
		$gpid = 82;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 4/m") {
		$gpid = 83;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 42/m") {
		$gpid = 84;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 4/n :1") {
		$gpid = 85;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 4/n :2") {
		$gpid = 85;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 42/n :1") {
		$gpid = 86;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 42/n :2") {
		$gpid = 86;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 4/m") {
		$gpid = 87;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 41/a :1") {
		$gpid = 88;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.75];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.25];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.75];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 41/a :2") {
		$gpid = 88;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.75,0.25,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.75,0.75,0.75];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.25,0.25,0.25];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 4 2 2") {
		$gpid = 89;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 4 21 2") {
		$gpid = 90;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 41 2 2") {
		$gpid = 91;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.75];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.25];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 41 21 2") {
		$gpid = 92;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.75];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.25];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 42 2 2") {
		$gpid = 93;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 42 21 2") {
		$gpid = 94;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 43 2 2") {
		$gpid = 95;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.75];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.25];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.25];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.75];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 43 21 2") {
		$gpid = 96;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.75];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.25];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.25];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.75];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 4 2 2") {
		$gpid = 97;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 41 2 2") {
		$gpid = 98;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.25];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.75];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 4 m m") {
		$gpid = 99;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 4 b m") {
		$gpid = 100;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 42 c m") {
		$gpid = 101;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 42 n m") {
		$gpid = 102;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 4 c c") {
		$gpid = 103;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 4 n c") {
		$gpid = 104;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 42 m c") {
		$gpid = 105;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P 42 b c") {
		$gpid = 106;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I 4 m m") {
		$gpid = 107;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I 4 c m") {
		$gpid = 108;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I 41 m d") {
		$gpid = 109;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.75];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.25];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.75];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "I 41 c d") {
		$gpid = 110;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.75];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.25];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,0] ];
	}
	elsif ($id eq "P -4 2 m") {
		$gpid = 111;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P -4 2 c") {
		$gpid = 112;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P -4 21 m") {
		$gpid = 113;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P -4 21 c") {
		$gpid = 114;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P -4 m 2") {
		$gpid = 115;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P -4 c 2") {
		$gpid = 116;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P -4 b 2") {
		$gpid = 117;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P -4 n 2") {
		$gpid = 118;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I -4 m 2") {
		$gpid = 119;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I -4 c 2") {
		$gpid = 120;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I -4 2 m") {
		$gpid = 121;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I -4 2 d") {
		$gpid = 122;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.25];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.75];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.25];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.75];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 4/m m m") {
		$gpid = 123;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 4/m c c") {
		$gpid = 124;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 4/n b m :1") {
		$gpid = 125;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 4/n b m :2") {
		$gpid = 125;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 4/n n c :1") {
		$gpid = 126;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 4/n n c :2") {
		$gpid = 126;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 4/m b m") {
		$gpid = 127;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 4/m n c") {
		$gpid = 128;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 4/n m m :1") {
		$gpid = 129;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 4/n m m :2") {
		$gpid = 129;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 4/n c c :1") {
		$gpid = 130;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 4/n c c :2") {
		$gpid = 130;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 42/m m c") {
		$gpid = 131;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 42/m c m") {
		$gpid = 132;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 42/n b c :1") {
		$gpid = 133;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 42/n b c :2") {
		$gpid = 133;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 42/n n m :1") {
		$gpid = 134;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 42/n n m :2") {
		$gpid = 134;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 42/m b c") {
		$gpid = 135;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 42/m n m") {
		$gpid = 136;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 42/n m c :1") {
		$gpid = 137;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 42/n m c :2") {
		$gpid = 137;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 42/n c m :1") {
		$gpid = 138;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 42/n c m :2") {
		$gpid = 138;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 4/m m m") {
		$gpid = 139;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 4/m c m") {
		$gpid = 140;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 41/a m d :1") {
		$gpid = 141;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.25];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.75];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.25];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.75];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.25];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 41/a m d :2") {
		$gpid = 141;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.25,0.75,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.25,0.25,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.25,0.75,0.25];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.25,0.25,0.75];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.75,0.75,0.25];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.75,0.75,0.25];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 41/a c d :1") {
		$gpid = 142;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.75];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.25];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.25];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.75];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.25];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.75];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 41/a c d :2") {
		$gpid = 142;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.25,0.75,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.25,0.25,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.75,0.75,0.25];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.75,0.25,0.25];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.75,0.75,0.75];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "P 3") {
		$gpid = 143;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,0] ];
	}
	elsif ($id eq "P 31") {
		$gpid = 144;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Cs, [0,0,0];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,0] ];
	}
	elsif ($id eq "P 32") {
		$gpid = 145;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Cs, [0,0,0];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,0] ];
	}
	elsif ($id eq "R 3 :H") {
		$gpid = 146;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.666666666666667,0.333333333333333,0.333333333333333];
		push @Cs, [0.333333333333333,0.666666666666667,0.666666666666667];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,0] ];
	}
	elsif ($id eq "R 3 :R") {
		$gpid = 146;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,0] ];
	}
	elsif ($id eq "P -3") {
		$gpid = 147;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "R -3 :H") {
		$gpid = 148;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.666666666666667,0.333333333333333,0.333333333333333];
		push @Cs, [0.333333333333333,0.666666666666667,0.666666666666667];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "R -3 :R") {
		$gpid = 148;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "P 3 1 2") {
		$gpid = 149;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,1/2] ];
	}
	elsif ($id eq "P 3 2 1") {
		$gpid = 150;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "P 31 1 2") {
		$gpid = 151;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Cs, [0,0,0];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,1/2] ];
	}
	elsif ($id eq "P 31 2 1") {
		$gpid = 152;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "P 32 1 2") {
		$gpid = 153;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Cs, [0,0,0];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,1/2] ];
	}
	elsif ($id eq "P 32 2 1") {
		$gpid = 154;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "R 3 2 :H") {
		$gpid = 155;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.666666666666667,0.333333333333333,0.333333333333333];
		push @Cs, [0.333333333333333,0.666666666666667,0.666666666666667];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "R 3 2 :R") {
		$gpid = 155;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "P 3 m 1") {
		$gpid = 156;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,0] ];
	}
	elsif ($id eq "P 3 1 m") {
		$gpid = 157;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,0] ];
	}
	elsif ($id eq "P 3 c 1") {
		$gpid = 158;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,0] ];
	}
	elsif ($id eq "P 3 1 c") {
		$gpid = 159;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,0] ];
	}
	elsif ($id eq "R 3 m :H") {
		$gpid = 160;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.666666666666667,0.333333333333333,0.333333333333333];
		push @Cs, [0.333333333333333,0.666666666666667,0.666666666666667];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,0] ];
	}
	elsif ($id eq "R 3 m :R") {
		$gpid = 160;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,0] ];
	}
	elsif ($id eq "R 3 c :H") {
		$gpid = 161;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.666666666666667,0.333333333333333,0.333333333333333];
		push @Cs, [0.333333333333333,0.666666666666667,0.666666666666667];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,0] ];
	}
	elsif ($id eq "R 3 c :R") {
		$gpid = 161;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,0] ];
	}
	elsif ($id eq "P -3 1 m") {
		$gpid = 162;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "P -3 1 c") {
		$gpid = 163;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "P -3 m 1") {
		$gpid = 164;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "P -3 c 1") {
		$gpid = 165;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "R -3 m :H") {
		$gpid = 166;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.666666666666667,0.333333333333333,0.333333333333333];
		push @Cs, [0.333333333333333,0.666666666666667,0.666666666666667];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "R -3 m :R") {
		$gpid = 166;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "R -3 c :H") {
		$gpid = 167;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.666666666666667,0.333333333333333,0.333333333333333];
		push @Cs, [0.333333333333333,0.666666666666667,0.666666666666667];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "R -3 c :R") {
		$gpid = 167;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "P 6") {
		$gpid = 168;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,0] ];
	}
	elsif ($id eq "P 61") {
		$gpid = 169;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.166666666666667];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.833333333333333];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,0] ];
	}
	elsif ($id eq "P 65") {
		$gpid = 170;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.833333333333333];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.166666666666667];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,0] ];
	}
	elsif ($id eq "P 62") {
		$gpid = 171;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,0] ];
	}
	elsif ($id eq "P 64") {
		$gpid = 172;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,0] ];
	}
	elsif ($id eq "P 63") {
		$gpid = 173;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,0] ];
	}
	elsif ($id eq "P -6") {
		$gpid = 174;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,1/2] ];
	}
	elsif ($id eq "P 6/m") {
		$gpid = 175;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "P 63/m") {
		$gpid = 176;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "P 6 2 2") {
		$gpid = 177;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "P 61 2 2") {
		$gpid = 178;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.166666666666667];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.833333333333333];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.833333333333333];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.166666666666667];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "P 65 2 2") {
		$gpid = 179;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.833333333333333];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.166666666666667];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.166666666666667];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.833333333333333];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "P 62 2 2") {
		$gpid = 180;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "P 64 2 2") {
		$gpid = 181;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.333333333333333];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.666666666666667];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "P 63 2 2") {
		$gpid = 182;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1/2] ];
	}
	elsif ($id eq "P 6 m m") {
		$gpid = 183;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,0] ];
	}
	elsif ($id eq "P 6 c c") {
		$gpid = 184;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,0] ];
	}
	elsif ($id eq "P 63 c m") {
		$gpid = 185;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,0] ];
	}
	elsif ($id eq "P 63 m c") {
		$gpid = 186;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,0] ];
	}
	elsif ($id eq "P -6 m 2") {
		$gpid = 187;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,1/2] ];
	}
	elsif ($id eq "P -6 c 2") {
		$gpid = 188;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,1/2] ];
	}
	elsif ($id eq "P -6 2 m") {
		$gpid = 189;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,1/2] ];
	}
	elsif ($id eq "P -6 2 c") {
		$gpid = 190;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,1/2] ];
	}
	elsif ($id eq "P 6/m m m") {
		$gpid = 191;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,1/2] ];
	}
	elsif ($id eq "P 6/m c c") {
		$gpid = 192;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,1/2] ];
	}
	elsif ($id eq "P 63/m c m") {
		$gpid = 193;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,1/2] ];
	}
	elsif ($id eq "P 63/m m c") {
		$gpid = 194;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [-1,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,1,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [-1,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,-1,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [1,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,2/3] , [0,2/3] , [0,1/2] ];
	}
	elsif ($id eq "P 2 3") {
		$gpid = 195;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "F 2 3") {
		$gpid = 196;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 2 3") {
		$gpid = 197;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "P 21 3") {
		$gpid = 198;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "I 21 3") {
		$gpid = 199;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "P m -3") {
		$gpid = 200;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "P n -3 :1") {
		$gpid = 201;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "P n -3 :2") {
		$gpid = 201;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "F m -3") {
		$gpid = 202;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "F d -3 :1") {
		$gpid = 203;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.75,0.75,0.25];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.75,0.75,0.25];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.75,0.75,0.25];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "F d -3 :2") {
		$gpid = 203;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.25,0.25,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.25,0.25];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.25,0,0.25];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.25,0.25,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0.25,0.25];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.25,0,0.25];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0.25,0.25];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.25,0,0.25];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.25,0.25,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.75,0.75,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.75,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.75,0,0.75];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.75,0.75,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0.75,0.75];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.75,0,0.75];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0.75,0.75];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.75,0,0.75];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.75,0.75,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I m -3") {
		$gpid = 204;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "P a -3") {
		$gpid = 205;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "I a -3") {
		$gpid = 206;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "P 4 3 2") {
		$gpid = 207;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "P 42 3 2") {
		$gpid = 208;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "F 4 3 2") {
		$gpid = 209;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "F 41 3 2") {
		$gpid = 210;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.75,0.75,0.25];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I 4 3 2") {
		$gpid = 211;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "P 43 3 2") {
		$gpid = 212;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.75,0.75,0.25];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.75,0.75,0.25];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.75,0.75,0.25];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.75,0.25,0.75];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "P 41 3 2") {
		$gpid = 213;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.25,0.75,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.25,0.25,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.75,0.25,0.25];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.75,0.75,0.75];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.25,0.75,0.25];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.25,0.25,0.75];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.75,0.25,0.25];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.75,0.75,0.75];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.75,0.25,0.25];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.75,0.75,0.75];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.25,0.25,0.75];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.25,0.75,0.25];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "I 41 3 2") {
		$gpid = 214;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.25,0.75,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.25,0.25,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.25,0.75,0.25];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.25,0.25,0.75];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.75,0.25,0.25];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.75,0.75,0.25];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.75,0.25,0.75];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "P -4 3 m") {
		$gpid = 215;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "F -4 3 m") {
		$gpid = 216;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I -4 3 m") {
		$gpid = 217;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "P -4 3 n") {
		$gpid = 218;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "F -4 3 c") {
		$gpid = 219;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,0,1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I -4 3 d") {
		$gpid = 220;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.25,0.75,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.75,0.75,0.25];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.75,0.75,0.75];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.25,0.75,0.25];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.75,0.75,0.25];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.75,0.75,0.75];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,-1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.25,0.25,0.75];
		push @Rs, [ [0,0,1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.75,0.25,0.75];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "P m -3 m") {
		$gpid = 221;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "P n -3 n :1") {
		$gpid = 222;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "P n -3 n :2") {
		$gpid = 222;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "P m -3 n") {
		$gpid = 223;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "P n -3 m :1") {
		$gpid = 224;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,-1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "P n -3 m :2") {
		$gpid = 224;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Cs, [0,0,0];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "F m -3 m") {
		$gpid = 225;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "F m -3 c") {
		$gpid = 226;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,0,1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "F d -3 m :1") {
		$gpid = 227;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.75,0.75,0.25];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.75,0.75,0.25];
		push @Rs, [ [0,0,-1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "F d -3 m :2") {
		$gpid = 227;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0.25,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.75,0.25,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.75,0,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.25,0.25];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.75,0.25,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.75,0,0.75];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0.25,0.25];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.75,0.25,0.5];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.75,0,0.75];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0.25,0.25];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.75,0.25,0.5];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.75,0,0.75];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0.5,0.75,0.25];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.25,0.75,0.5];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.25,0.5,0.75];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.25,0.25,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.25,0,0.25];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.5,0.25,0.75];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0.75,0.75];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.25,0.75,0.5];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.25,0,0.25];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.75,0.75];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.25,0.75,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.25,0,0.25];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0.75,0.75];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.25,0.75,0.5];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.25,0,0.25];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0.75,0.75];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.25,0.75,0.5];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.25,0,0.25];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0.5,0.25,0.75];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.75,0.25,0.5];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.75,0.5,0.25];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.75,0.75,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.75,0,0.75];
		push @Rs, [ [0,0,1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.5,0.75,0.25];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "F d -3 c :1") {
		$gpid = 228;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.75,0.75,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0.75,0.75,0.75];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.75,0.25,0.25];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0.75,0.75,0.75];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.75,0.25,0.25];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0.75,0.75,0.75];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.75,0.25,0.25];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0.75,0.75,0.75];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.75,0.25,0.25];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0.75,0.75,0.75];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0.25,0.75,0.25];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.25,0.25,0.75];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.75,0.75,0.75];
		push @Rs, [ [0,0,-1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "F d -3 c :2") {
		$gpid = 228;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.5,0.25,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.25,0.75,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.75,0.5,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0.25,0.25];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.25,0.25,0.5];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.25,0.5,0.75];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.5,0.25,0.25];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.25,0.75,0.5];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.75,0.5,0.75];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0.25,0.25];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.25,0.25,0.5];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.25,0.5,0.75];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0.5,0.25,0.75];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.25,0.75,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.75,0.5,0.25];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.25,0.25,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.25,0,0.75];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0.75,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.5,0.75,0.75];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.75,0.25,0.5];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.25,0.5,0.25];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0.75,0.75];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.75,0.75,0.5];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.75,0.5,0.25];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.5,0.75,0.75];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.75,0.25,0.5];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.25,0.5,0.25];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0.75,0.75];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.75,0.75,0.5];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.75,0.5,0.25];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0.5,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0.5,0.75,0.25];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.75,0.25,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.25,0.5,0.75];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.5,0.5,0.5];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.75,0.75,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.75,0,0.25];
		push @Rs, [ [0,0,1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0.25,0.75];
		push @Cs, [0,0,0];
		push @Cs, [0,0.5,0.5];
		push @Cs, [0.5,0,0.5];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1/2] , [0,1/2] , [0,1/2] ];
	}
	elsif ($id eq "I m -3 m") {
		$gpid = 229;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,-1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,0,1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "I a -3 d") {
		$gpid = 230;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.25,0.75,0.25];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.25,0.25,0.75];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.25,0.75,0.25];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.25,0.25,0.75];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.75,0.25,0.25];
		push @Rs, [ [0,-1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.25,0.25,0.25];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.75,0.75,0.25];
		push @Rs, [ [0,0,-1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [-1,0,0] , [0,-1,0] , [0,0,-1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,1,0] , [-1,0,0] , [0,0,-1] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,-1,0] , [1,0,0] , [0,0,-1] ];
		push @Ts, [0.75,0.75,0.25];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [0,-1,0] , [-1,0,0] , [0,0,1] ];
		push @Ts, [0.75,0.25,0.25];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [0,1,0] , [1,0,0] , [0,0,1] ];
		push @Ts, [0.75,0.75,0.75];
		push @Rs, [ [0,0,-1] , [-1,0,0] , [0,-1,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,0,-1] , [0,-1,0] ];
		push @Ts, [0.75,0.25,0.75];
		push @Rs, [ [0,0,1] , [1,0,0] , [0,-1,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [-1,0,0] , [0,0,1] , [0,-1,0] ];
		push @Ts, [0.75,0.75,0.25];
		push @Rs, [ [0,0,-1] , [1,0,0] , [0,1,0] ];
		push @Ts, [0,0,0.5];
		push @Rs, [ [-1,0,0] , [0,0,-1] , [0,1,0] ];
		push @Ts, [0.75,0.25,0.25];
		push @Rs, [ [0,0,1] , [-1,0,0] , [0,1,0] ];
		push @Ts, [0.5,0,0];
		push @Rs, [ [1,0,0] , [0,0,1] , [0,1,0] ];
		push @Ts, [0.75,0.75,0.75];
		push @Rs, [ [0,-1,0] , [0,0,-1] , [-1,0,0] ];
		push @Ts, [0,0,0];
		push @Rs, [ [0,-1,0] , [0,0,1] , [1,0,0] ];
		push @Ts, [0.5,0.5,0];
		push @Rs, [ [0,0,-1] , [0,-1,0] , [1,0,0] ];
		push @Ts, [0.25,0.75,0.75];
		push @Rs, [ [0,1,0] , [0,0,-1] , [1,0,0] ];
		push @Ts, [0,0.5,0.5];
		push @Rs, [ [0,0,1] , [0,1,0] , [1,0,0] ];
		push @Ts, [0.75,0.75,0.75];
		push @Rs, [ [0,1,0] , [0,0,1] , [-1,0,0] ];
		push @Ts, [0.5,0,0.5];
		push @Rs, [ [0,0,-1] , [0,1,0] , [-1,0,0] ];
		push @Ts, [0.25,0.25,0.75];
		push @Rs, [ [0,0,1] , [0,-1,0] , [-1,0,0] ];
		push @Ts, [0.25,0.75,0.25];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "I 1 21 1") {
		$gpid = 5;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0.5];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0.5];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "C 1 21 1") {
		$gpid = 5;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [-1,0,0] , [0,1,0] , [0,0,-1] ];
		push @Ts, [0.5,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	elsif ($id eq "B 1 1 m") {
		$gpid = 8;
		push @Rs, [ [1,0,0] , [0,1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Rs, [ [1,0,0] , [0,-1,0] , [0,0,1] ];
		push @Ts, [0,0,0];
		push @Cs, [0,0,0];
		push @Cs, [0.5,0.5,0];
		$cheshire = [ [0,1] , [0,1] , [0,1] ];
	}
	else {
		print STDERR "SPACEGROUP $id undefined!\n";
		exit -1;
	}

	# expand cenops
	foreach my $C (@Cs) {
		foreach my $i (0..$#Ts) {
			push @R_all, $Rs[$i];
			push @T_all, vadd( $Ts[ $i ] , $C );
		}
	}
	my $nsymm = $#T_all + 1;

	return ($gpid, $nsymm, \@R_all, \@T_all, $cheshire);
}

#######
# end #
#######

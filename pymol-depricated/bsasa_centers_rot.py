import pymol
from pymol import cmd
import sys,os,random,re
from pymol.cgo import *
import random

POCKET1 = """
@subgroup {    Pocket1} dominant off
@vectorlist {Pck1} color=green
{ ca  arg      7 } P    2.803,   8.648,   4.367 { ca  arg      7 } L    2.803,   8.648,   4.367
{ cd1 leu     10 } L    1.867,   7.291,   9.151
{ cb  lys     48 } P   -6.083,   8.509,  14.663 { cb  lys     48 } L   -6.083,   8.509,  14.663
{ cg  leu     49 } L   -1.484,  10.094,  15.959
{ o   asp     45 } P   -4.102,  11.719,  14.877 { o   asp     45 } L   -4.102,  11.719,  14.877
{ cg  leu     49 } L   -1.484,  10.094,  15.959
{ o   asp     45 } P   -4.102,  11.719,  14.877 { o   asp     45 } L   -4.102,  11.719,  14.877
{ cb  lys     48 } L   -6.083,   8.509,  14.663
{ cg  leu     49 } P   -1.484,  10.094,  15.959 { cg  leu     49 } L   -1.484,  10.094,  15.959
{ cd1 leu     49 } L   -0.561,  11.234,  15.551
{ cb  asp     45 } P   -3.836,  13.216,  11.899 { cb  asp     45 } L   -3.836,  13.216,  11.899
{ cd1 leu     49 } L   -0.561,  11.234,  15.551
{ cb  asp     45 } P   -3.836,  13.216,  11.899 { cb  asp     45 } L   -3.836,  13.216,  11.899
{ cg  leu     49 } L   -1.484,  10.094,  15.959
{ ca  asp     45 } P   -4.884,  12.484,  12.742 { ca  asp     45 } L   -4.884,  12.484,  12.742
{ cg  leu     49 } L   -1.484,  10.094,  15.959
{ ca  asp     45 } P   -4.884,  12.484,  12.742 { ca  asp     45 } L   -4.884,  12.484,  12.742
{ cb  asp     45 } L   -3.836,  13.216,  11.899
{ od1 asp     45 } P   -4.455,  11.853,  10.086 { od1 asp     45 } L   -4.455,  11.853,  10.086
{ od2 asp     45 } L   -3.344,  13.600,   9.633
{ ce2 phe     37 } P   -3.372,  14.257,   5.145 { ce2 phe     37 } L   -3.372,  14.257,   5.145
{ od2 asp     45 } L   -3.344,  13.600,   9.633
{ ce2 phe     37 } P   -3.372,  14.257,   5.145 { ce2 phe     37 } L   -3.372,  14.257,   5.145
{ od1 asp     45 } L   -4.455,  11.853,  10.086
{ cd1 leu     10 } P    1.867,   7.291,   9.151 { cd1 leu     10 } L    1.867,   7.291,   9.151
{ cd2 leu     10 } L    1.728,   8.986,  10.991
{ cg2 val      9 } P    0.064,  14.243,   6.530 { cg2 val      9 } L    0.064,  14.243,   6.530
{ ce2 phe     37 } L   -3.372,  14.257,   5.145
{ ca  gln      6 } P   -0.815,   9.814,   4.283 { ca  gln      6 } L   -0.815,   9.814,   4.283
{ ce2 phe     37 } L   -3.372,  14.257,   5.145
{ ca  gln      6 } P   -0.815,   9.814,   4.283 { ca  gln      6 } L   -0.815,   9.814,   4.283
{ cg2 val      9 } L    0.064,  14.243,   6.530
{ cd1 leu     10 } P    1.867,   7.291,   9.151 { cd1 leu     10 } L    1.867,   7.291,   9.151
{ od1 asp     45 } L   -4.455,  11.853,  10.086
{ cg  leu     10 } P    2.055,   8.756,   9.523 { cg  leu     10 } L    2.055,   8.756,   9.523
{ od1 asp     45 } L   -4.455,  11.853,  10.086
{ cg  leu     10 } P    2.055,   8.756,   9.523 { cg  leu     10 } L    2.055,   8.756,   9.523
{ cd1 leu     10 } L    1.867,   7.291,   9.151
{ cb  gln      6 } P   -1.660,   8.539,   4.337 { cb  gln      6 } L   -1.660,   8.539,   4.337
{ od1 asp     45 } L   -4.455,  11.853,  10.086
{ cb  gln      6 } P   -1.660,   8.539,   4.337 { cb  gln      6 } L   -1.660,   8.539,   4.337
{ cd1 leu     10 } L    1.867,   7.291,   9.151
{ cb  gln      6 } P   -1.660,   8.539,   4.337 { cb  gln      6 } L   -1.660,   8.539,   4.337
{ cg  leu     10 } L    2.055,   8.756,   9.523
{ cg  asp     45 } P   -3.833,  12.831,  10.426 { cg  asp     45 } L   -3.833,  12.831,  10.426
{ od1 asp     45 } L   -4.455,  11.853,  10.086
{ cg1 val      9 } P    0.536,  13.866,   8.955 { cg1 val      9 } L    0.536,  13.866,   8.955
{ od1 asp     45 } L   -4.455,  11.853,  10.086
{ cg1 val      9 } P    0.536,  13.866,   8.955 { cg1 val      9 } L    0.536,  13.866,   8.955
{ cg  asp     45 } L   -3.833,  12.831,  10.426
{ cb  val      9 } P    0.829,  13.398,   7.537 { cb  val      9 } L    0.829,  13.398,   7.537
{ od1 asp     45 } L   -4.455,  11.853,  10.086
{ cb  val      9 } P    0.829,  13.398,   7.537 { cb  val      9 } L    0.829,  13.398,   7.537
{ cg  asp     45 } L   -3.833,  12.831,  10.426
{ cb  val      9 } P    0.829,  13.398,   7.537 { cb  val      9 } L    0.829,  13.398,   7.537
{ cg1 val      9 } L    0.536,  13.866,   8.955
{ cg2 val      9 } P    0.064,  14.243,   6.530 { cg2 val      9 } L    0.064,  14.243,   6.530
{ od2 asp     45 } L   -3.344,  13.600,   9.633
{ cb  val      9 } P    0.829,  13.398,   7.537 { cb  val      9 } L    0.829,  13.398,   7.537
{ od2 asp     45 } L   -3.344,  13.600,   9.633
{ cb  val      9 } P    0.829,  13.398,   7.537 { cb  val      9 } L    0.829,  13.398,   7.537
{ ce2 phe     37 } L   -3.372,  14.257,   5.145
{ cb  val      9 } P    0.829,  13.398,   7.537 { cb  val      9 } L    0.829,  13.398,   7.537
{ cg2 val      9 } L    0.064,  14.243,   6.530
{ ce  lys     48 } P   -4.977,   5.372,  12.691 { ce  lys     48 } L   -4.977,   5.372,  12.691
{ cd2 leu     49 } L   -0.765,   8.757,  15.863
{ cd  lys     48 } P   -5.427,   6.804,  12.936 { cd  lys     48 } L   -5.427,   6.804,  12.936
{ cd2 leu     49 } L   -0.765,   8.757,  15.863
{ cd  lys     48 } P   -5.427,   6.804,  12.936 { cd  lys     48 } L   -5.427,   6.804,  12.936
{ ce  lys     48 } L   -4.977,   5.372,  12.691
{ cd2 leu     10 } P    1.728,   8.986,  10.991 { cd2 leu     10 } L    1.728,   8.986,  10.991
{ cd2 leu     49 } L   -0.765,   8.757,  15.863
{ cd2 leu     10 } P    1.728,   8.986,  10.991 { cd2 leu     10 } L    1.728,   8.986,  10.991
{ ce  lys     48 } L   -4.977,   5.372,  12.691
{ cd2 leu     10 } P    1.728,   8.986,  10.991 { cd2 leu     10 } L    1.728,   8.986,  10.991
{ cd  lys     48 } L   -5.427,   6.804,  12.936
{ cd2 leu     49 } P   -0.765,   8.757,  15.863 { cd2 leu     49 } L   -0.765,   8.757,  15.863
{ ce2 tyr     74 } L    3.275,   7.836,  14.253
{ cd1 leu     52 } P   -3.221,   4.984,  16.696 { cd1 leu     52 } L   -3.221,   4.984,  16.696
{ ce2 phe     73 } L    1.133,   3.773,  14.096
{ cd2 leu     49 } P   -0.765,   8.757,  15.863 { cd2 leu     49 } L   -0.765,   8.757,  15.863
{ ce2 phe     73 } L    1.133,   3.773,  14.096
{ cd2 leu     49 } P   -0.765,   8.757,  15.863 { cd2 leu     49 } L   -0.765,   8.757,  15.863
{ cd1 leu     52 } L   -3.221,   4.984,  16.696
{ ce  lys     48 } P   -4.977,   5.372,  12.691 { ce  lys     48 } L   -4.977,   5.372,  12.691
{ ce2 phe     73 } L    1.133,   3.773,  14.096
{ ce  lys     48 } P   -4.977,   5.372,  12.691 { ce  lys     48 } L   -4.977,   5.372,  12.691
{ cd1 leu     52 } L   -3.221,   4.984,  16.696
{ cg1 val      9 } P    0.536,  13.866,   8.955 { cg1 val      9 } L    0.536,  13.866,   8.955
{ cd1 ile     13 } L    2.644,  12.675,  12.899
{ o   gln      6 } P    0.873,   9.804,   5.988 { o   gln      6 } L    0.873,   9.804,   5.988
{ cd1 leu     10 } L    1.867,   7.291,   9.151
{ c   gln      6 } P    0.580,   9.529,   4.825 { c   gln      6 } L    0.580,   9.529,   4.825
{ cd1 leu     10 } L    1.867,   7.291,   9.151
{ cb  gln      6 } P   -1.660,   8.539,   4.337 { cb  gln      6 } L   -1.660,   8.539,   4.337
{ o   gln      6 } L    0.873,   9.804,   5.988
{ cb  gln      6 } P   -1.660,   8.539,   4.337 { cb  gln      6 } L   -1.660,   8.539,   4.337
{ c   gln      6 } L    0.580,   9.529,   4.825
{ cd  gln      6 } P   -3.929,   7.466,   4.075 { cd  gln      6 } L   -3.929,   7.466,   4.075
{ od1 asp     45 } L   -4.455,  11.853,  10.086
{ cg  gln      6 } P   -3.140,   8.762,   4.076 { cg  gln      6 } L   -3.140,   8.762,   4.076
{ od1 asp     45 } L   -4.455,  11.853,  10.086
{ cg  gln      6 } P   -3.140,   8.762,   4.076 { cg  gln      6 } L   -3.140,   8.762,   4.076
{ cd  gln      6 } L   -3.929,   7.466,   4.075
{ cd1 ile     13 } P    2.644,  12.675,  12.899 { cd1 ile     13 } L    2.644,  12.675,  12.899
{ cb  asp     45 } L   -3.836,  13.216,  11.899
{ cb  gln      6 } P   -1.660,   8.539,   4.337 { cb  gln      6 } L   -1.660,   8.539,   4.337
{ oe1 gln      6 } L   -3.369,   6.382,   4.261
{ cd2 leu      2 } P   -1.227,   2.531,   4.222 { cd2 leu      2 } L   -1.227,   2.531,   4.222
{ oe1 gln      6 } L   -3.369,   6.382,   4.261
{ ca  gln      6 } P   -0.815,   9.814,   4.283 { ca  gln      6 } L   -0.815,   9.814,   4.283
{ cg  gln      6 } L   -3.140,   8.762,   4.076
{ cg  leu     10 } P    2.055,   8.756,   9.523 { cg  leu     10 } L    2.055,   8.756,   9.523
{ cd2 leu     10 } L    1.728,   8.986,  10.991
{ cd1 leu     10 } P    1.867,   7.291,   9.151 { cd1 leu     10 } L    1.867,   7.291,   9.151
{ nz  lys     48 } L   -4.743,   5.101,  11.247
{ oe1 gln      6 } P   -3.369,   6.382,   4.261 { oe1 gln      6 } L   -3.369,   6.382,   4.261
{ nz  lys     48 } L   -4.743,   5.101,  11.247
{ oe1 gln      6 } P   -3.369,   6.382,   4.261 { oe1 gln      6 } L   -3.369,   6.382,   4.261
{ cd1 leu     10 } L    1.867,   7.291,   9.151
{ cb  gln      6 } P   -1.660,   8.539,   4.337 { cb  gln      6 } L   -1.660,   8.539,   4.337
{ nz  lys     48 } L   -4.743,   5.101,  11.247
{ od1 asp     45 } P   -4.455,  11.853,  10.086 { od1 asp     45 } L   -4.455,  11.853,  10.086
{ cd  lys     48 } L   -5.427,   6.804,  12.936
{ cd  gln      6 } P   -3.929,   7.466,   4.075 { cd  gln      6 } L   -3.929,   7.466,   4.075
{ nz  lys     48 } L   -4.743,   5.101,  11.247
{ cd  gln      6 } P   -3.929,   7.466,   4.075 { cd  gln      6 } L   -3.929,   7.466,   4.075
{ oe1 gln      6 } L   -3.369,   6.382,   4.261
{ cg  gln      6 } P   -3.140,   8.762,   4.076 { cg  gln      6 } L   -3.140,   8.762,   4.076
{ ce2 phe     37 } L   -3.372,  14.257,   5.145
{ c   gln      6 } P    0.580,   9.529,   4.825 { c   gln      6 } L    0.580,   9.529,   4.825
{ ca  arg      7 } L    2.803,   8.648,   4.367
{ cd1 leu      2 } P    1.057,   3.170,   5.027 { cd1 leu      2 } L    1.057,   3.170,   5.027
{ ca  arg      7 } L    2.803,   8.648,   4.367
{ cd1 leu      2 } P    1.057,   3.170,   5.027 { cd1 leu      2 } L    1.057,   3.170,   5.027
{ c   gln      6 } L    0.580,   9.529,   4.825
{ cd1 leu     49 } P   -0.561,  11.234,  15.551 { cd1 leu     49 } L   -0.561,  11.234,  15.551
{ cd2 leu     49 } L   -0.765,   8.757,  15.863
{ cd2 leu     10 } P    1.728,   8.986,  10.991 { cd2 leu     10 } L    1.728,   8.986,  10.991
{ ce2 tyr     74 } L    3.275,   7.836,  14.253
{ cd2 leu     10 } P    1.728,   8.986,  10.991 { cd2 leu     10 } L    1.728,   8.986,  10.991
{ cd1 leu     49 } L   -0.561,  11.234,  15.551
{ nz  lys     48 } P   -4.743,   5.101,  11.247 { nz  lys     48 } L   -4.743,   5.101,  11.247
{ ce2 phe     73 } L    1.133,   3.773,  14.096
{ ce  lys     48 } P   -4.977,   5.372,  12.691 { ce  lys     48 } L   -4.977,   5.372,  12.691
{ nz  lys     48 } L   -4.743,   5.101,  11.247
{ cd1 leu     10 } P    1.867,   7.291,   9.151 { cd1 leu     10 } L    1.867,   7.291,   9.151
{ ce  lys     48 } L   -4.977,   5.372,  12.691
{ cg  gln      6 } P   -3.140,   8.762,   4.076 { cg  gln      6 } L   -3.140,   8.762,   4.076
{ nz  lys     48 } L   -4.743,   5.101,  11.247
{ cb  gln      6 } P   -1.660,   8.539,   4.337 { cb  gln      6 } L   -1.660,   8.539,   4.337
{ cd  gln      6 } L   -3.929,   7.466,   4.075
{ cb  gln      6 } P   -1.660,   8.539,   4.337 { cb  gln      6 } L   -1.660,   8.539,   4.337
{ cg  gln      6 } L   -3.140,   8.762,   4.076
{ cd2 leu     10 } P    1.728,   8.986,  10.991 { cd2 leu     10 } L    1.728,   8.986,  10.991
{ cg  leu     49 } L   -1.484,  10.094,  15.959
{ cd2 leu     10 } P    1.728,   8.986,  10.991 { cd2 leu     10 } L    1.728,   8.986,  10.991
{ cb  asp     45 } L   -3.836,  13.216,  11.899
{ cg1 val      9 } P    0.536,  13.866,   8.955 { cg1 val      9 } L    0.536,  13.866,   8.955
{ od2 asp     45 } L   -3.344,  13.600,   9.633
{ cg1 val      9 } P    0.536,  13.866,   8.955 { cg1 val      9 } L    0.536,  13.866,   8.955
{ cb  asp     45 } L   -3.836,  13.216,  11.899
{ cb  lys     48 } P   -6.083,   8.509,  14.663 { cb  lys     48 } L   -6.083,   8.509,  14.663
{ cd  lys     48 } L   -5.427,   6.804,  12.936
{ o   gln      6 } P    0.873,   9.804,   5.988 { o   gln      6 } L    0.873,   9.804,   5.988
{ cb  val      9 } L    0.829,  13.398,   7.537
{ ca  asp     45 } P   -4.884,  12.484,  12.742 { ca  asp     45 } L   -4.884,  12.484,  12.742
{ od1 asp     45 } L   -4.455,  11.853,  10.086
{ cg  leu     10 } P    2.055,   8.756,   9.523 { cg  leu     10 } L    2.055,   8.756,   9.523
{ cg  asp     45 } L   -3.833,  12.831,  10.426
{ cg1 val      9 } P    0.536,  13.866,   8.955 { cg1 val      9 } L    0.536,  13.866,   8.955
{ cg  leu     10 } L    2.055,   8.756,   9.523
{ cd2 leu     10 } P    1.728,   8.986,  10.991 { cd2 leu     10 } L    1.728,   8.986,  10.991
{ cd1 ile     13 } L    2.644,  12.675,  12.899
{ cd  lys     48 } P   -5.427,   6.804,  12.936 { cd  lys     48 } L   -5.427,   6.804,  12.936
{ cg  leu     49 } L   -1.484,  10.094,  15.959
{ cb  asp     45 } P   -3.836,  13.216,  11.899 { cb  asp     45 } L   -3.836,  13.216,  11.899
{ cd  lys     48 } L   -5.427,   6.804,  12.936
{ ca  asp     45 } P   -4.884,  12.484,  12.742 { ca  asp     45 } L   -4.884,  12.484,  12.742
{ cd  lys     48 } L   -5.427,   6.804,  12.936
{ cb  val      9 } P    0.829,  13.398,   7.537 { cb  val      9 } L    0.829,  13.398,   7.537
{ cg  leu     10 } L    2.055,   8.756,   9.523
{ cd1 leu      2 } P    1.057,   3.170,   5.027 { cd1 leu      2 } L    1.057,   3.170,   5.027
{ nz  lys     48 } L   -4.743,   5.101,  11.247
{ cg  leu     49 } P   -1.484,  10.094,  15.959 { cg  leu     49 } L   -1.484,  10.094,  15.959
{ cd2 leu     49 } L   -0.765,   8.757,  15.863
{ ca  asp     45 } P   -4.884,  12.484,  12.742 { ca  asp     45 } L   -4.884,  12.484,  12.742
{ cb  lys     48 } L   -6.083,   8.509,  14.663
{ cd2 leu      2 } P   -1.227,   2.531,   4.222 { cd2 leu      2 } L   -1.227,   2.531,   4.222
{ nz  lys     48 } L   -4.743,   5.101,  11.247
{ cd1 leu      2 } P    1.057,   3.170,   5.027 { cd1 leu      2 } L    1.057,   3.170,   5.027
{ cd2 leu      2 } L   -1.227,   2.531,   4.222
{ ca  gln      6 } P   -0.815,   9.814,   4.283 { ca  gln      6 } L   -0.815,   9.814,   4.283
{ cb  gln      6 } L   -1.660,   8.539,   4.337
{ ce2 phe     73 } P    1.133,   3.773,  14.096 { ce2 phe     73 } L    1.133,   3.773,  14.096
{ ce2 tyr     74 } L    3.275,   7.836,  14.253
{ cd2 leu     10 } P    1.728,   8.986,  10.991 { cd2 leu     10 } L    1.728,   8.986,  10.991
{ ce2 phe     73 } L    1.133,   3.773,  14.096
{ o   gln      6 } P    0.873,   9.804,   5.988 { o   gln      6 } L    0.873,   9.804,   5.988
{ od1 asp     45 } L   -4.455,  11.853,  10.086
{ o   gln      6 } P    0.873,   9.804,   5.988 { o   gln      6 } L    0.873,   9.804,   5.988
{ cg  leu     10 } L    2.055,   8.756,   9.523
{ cg  asp     45 } P   -3.833,  12.831,  10.426 { cg  asp     45 } L   -3.833,  12.831,  10.426
{ od2 asp     45 } L   -3.344,  13.600,   9.633
{ ca  asp     45 } P   -4.884,  12.484,  12.742 { ca  asp     45 } L   -4.884,  12.484,  12.742
{ o   asp     45 } L   -4.102,  11.719,  14.877
{ ca  gln      6 } P   -0.815,   9.814,   4.283 { ca  gln      6 } L   -0.815,   9.814,   4.283
{ od1 asp     45 } L   -4.455,  11.853,  10.086
{ od1 asp     45 } P   -4.455,  11.853,  10.086 { od1 asp     45 } L   -4.455,  11.853,  10.086
{ nz  lys     48 } L   -4.743,   5.101,  11.247
{ ca  gln      6 } P   -0.815,   9.814,   4.283 { ca  gln      6 } L   -0.815,   9.814,   4.283
{ o   gln      6 } L    0.873,   9.804,   5.988
{ cd  lys     48 } P   -5.427,   6.804,  12.936 { cd  lys     48 } L   -5.427,   6.804,  12.936
{ nz  lys     48 } L   -4.743,   5.101,  11.247
{ cb  asp     45 } P   -3.836,  13.216,  11.899 { cb  asp     45 } L   -3.836,  13.216,  11.899
{ od1 asp     45 } L   -4.455,  11.853,  10.086
{ cb  asp     45 } P   -3.836,  13.216,  11.899 { cb  asp     45 } L   -3.836,  13.216,  11.899
{ cg  asp     45 } L   -3.833,  12.831,  10.426
{ ca  gln      6 } P   -0.815,   9.814,   4.283 { ca  gln      6 } L   -0.815,   9.814,   4.283
{ cb  val      9 } L    0.829,  13.398,   7.537
{ cd2 leu     10 } P    1.728,   8.986,  10.991 { cd2 leu     10 } L    1.728,   8.986,  10.991
{ cg  asp     45 } L   -3.833,  12.831,  10.426
{ cg1 val      9 } P    0.536,  13.866,   8.955 { cg1 val      9 } L    0.536,  13.866,   8.955
{ cd2 leu     10 } L    1.728,   8.986,  10.991
{ cd2 leu      2 } P   -1.227,   2.531,   4.222 { cd2 leu      2 } L   -1.227,   2.531,   4.222
{ cd1 leu     10 } L    1.867,   7.291,   9.151
{ cd1 leu      2 } P    1.057,   3.170,   5.027 { cd1 leu      2 } L    1.057,   3.170,   5.027
{ oe1 gln      6 } L   -3.369,   6.382,   4.261
{ cd1 leu      2 } P    1.057,   3.170,   5.027 { cd1 leu      2 } L    1.057,   3.170,   5.027
{ cb  gln      6 } L   -1.660,   8.539,   4.337
{ cd2 leu     10 } P    1.728,   8.986,  10.991 { cd2 leu     10 } L    1.728,   8.986,  10.991
{ od1 asp     45 } L   -4.455,  11.853,  10.086
1.867,   7.291,   9.151
"""

def loadcenters(file,nickname=None):
	if nickname is None:
		nickname = file.split('/')[-1][:4]
	print nickname
	cmd.load(file,nickname)
	cmd.do("useRosettaRadii")
	csel = nickname+"centers"
	psel = nickname+"rot"
	print nickname+" & resi 500-999"
	cmd.select(csel,nickname+" & resi 500-999")
	cmd.select(psel,nickname+" & resi   0-500")
	cmd.do("useOccRadii "+csel)
	cmd.hide('ev',nickname)
	cmd.show('sph',csel)
	cmd.show('cart',psel)
	cmd.show('line',psel)
	cmd.color('white',psel)

def loadpair(file):
	id = file.split('/')[-1][:4]
	loadcenters(file,id+'nat')
	loadcenters(file[:-17]+'_decoy_sasa_centers.pdb',id+'decoy')
	d = id+'decoy'
	n = id+'nat'
	cmd.align(n,d)

	

class LRvolume: 
	def __init__(self,line): # doesn't handle arbitrary order!
		rex = "^.*<volume"
		for tag in ("color","name","sa","v"):
			rex += "\s+%s=(\S+)"%tag
		rex += ".*>"		
		print rex
		self.color,self.name,self.sa,self.vol = re.match(rex,line).groups()

class LRarc:
	def __init__(self,line):
		rex = "^.*<arc"
		names = ("radius","theta_hi","theta_lo","x","y","z")
		for tag in names:
			rex += "\s+%s=(\S+)"%tag
		rex += "\s*>"		
		print rex
		vals = re.match(rex,line).groups()
		for i in range(len(names)):
			setattr(self,names[i],vals[i])
	def __repr__(self):
		return ",".join((self.radius,self.theta_hi,self.theta_lo,self.x,self.y,self.z))

s = "<volume color=white name=0 sa=4470.66 v=17019.7>"
a = "<arc radius=1.98063 theta_hi=180 theta_lo=0 x=-4.075 y=15.097 z=-6.52939> </arc>"

v = LRvolume(s)
r = LRarc(a)
print r

# { cb  val      9 } P    0.829,  13.398,   7.537 { cb  val      9 } L    0.829,  13.398,   7.537
# { cg1 val      9 } L    0.536,  13.866,   8.955

def parsekin(s):
	for l in s.split():
		pass
		

#parsekin(POCKET1)

#cmd.set("cgo_line_width",5)
#obj = [
#   BEGIN, LINES,
#   LINEWIDTH, 50,#
#
#   VERTEX, 2.803,   8.648,   4.367,
#   VERTEX, 1.867,   7.291,   9.151,#
#
#   END
#   ]                                                                                            
#
#cmd.load_cgo(obj,'cgo'+`random.random()`) 


def noclip():
	print "noclip!"
	cmd.clip("near",100)
	cmd.clip("far",-100)
	

cmd.extend("loadcenters",loadcenters)
cmd.extend("noclup",noclip)

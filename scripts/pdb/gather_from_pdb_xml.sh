#!/bin/bash
for fd in $@; do
	for f in $(find $fd -name \*.xml.gz); do
		code=$(zgrep      PDBx:pdbx_PDB_id_code                 $f | head -n1 | sed -e 's=^.*>\(.*\)<.*$=\1=g')
		resl=$(zgrep      PDBx:ls_d_res_high                    $f | head -n1 | sed -e 's=^.*>\(.*\)<.*$=\1=g')
		if [ ! $resl ]; then resl="0.00"; fi
		meth=$(zgrep "<PDBx:exptl entry_id=\"$code\" method=\"" $f | head -n1 | sed -e 's=^.*method\="\(.*\)">$=\1=g' | sed -e 's= =_=g')
		# printf "%4s %1.2f %10s\n" $code $resl $meth
		# echo $(echo CODE_$code | sed -e 's= =_=g') $(echo RESL_$resl | sed -e 's= =_=g') METH_$meth
		echo $(echo $code | sed -e 's= =_=g') $(echo $resl | sed -e 's= =_=g') $meth
	done
done

# resl ls_d_res_high resolution_high d_res_high
function runnewmatch {
    cd /work/sheffler/project/newmatch/fixdun
    for i in `cat scaff.list.rand$1`; do
	echo $i
	mkdir -p out/$i
	cd out/$i
	/usr/bin/time -f %E -o $i.time /work/sheffler/mini/bin/test_ikrs.omp.linuxiccrelease @input/ctp.flags -scaffold_active_site_residues ~/matchpos/${i}*.pos -s ~/matchpdb/${i}*.pdb.gz &> $i.log
	cd /work/sheffler/project/newmatch/fixdun
    done
}

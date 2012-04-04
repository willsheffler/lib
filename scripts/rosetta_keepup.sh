#!/bin/sh

function updaterosetta {
	cd $1/rosetta;
	svn up &> svn_up.log
	old=$(cat .revision)
	rev=$(svn info | grep 'Last Changed Rev:' | sed -e 's=Last Changed Rev: ==g')
	echo $rev > .revision
	if [[ $rev != $old ]]; then
		echo "change from rev $old to $rev"
		echo "cd $1/rosetta/rosetta_source"
		      cd $1/rosetta/rosetta_source
		echo "./update_options.sh &> autobuild/autobuild_opt_$rev.log"
		      ./update_options.sh &> autobuild/autobuild_opt_$rev.log
		echo "./scons.py --nover -j8 mode=release bin &> autobuild/autobuild_release_$rev.log"
		      ./scons.py --nover -j8 mode=release bin &> autobuild/autobuild_release_$rev.log
		echo "cd $1/rosetta/rosetta_tests/integration"
		      cd $1/rosetta/rosetta_tests/integration
		echo "python integration.py --database $1/rosetta/rosetta_database -j8 &> autobuild/integration_$rev.log"
		      python integration.py --database $1/rosetta/rosetta_database -j8 &> autobuild/integration_$rev.log
		echo "cd $1/rosetta/rosetta_source"
		      cd $1/rosetta/rosetta_source
		echo "./scons.py --nover -j8 mode=debug  &> autobuild/autobuild_debug_$rev.log"
		      ./scons.py --nover -j8 mode=debug  &> autobuild/autobuild_debug_$rev.log
		echo "./scons.py --nover -j8 cat=test &> autobuild/autobuild_test_$rev.log"
		      ./scons.py --nover -j8 cat=test &> autobuild/autobuild_test_$rev.log
		echo "python test/run.py --database $1/rosetta/rosetta_database -j8 &> autobuild/unittest_$rev.log"
		      python test/run.py --database $1/rosetta/rosetta_database -j8 &> autobuild/unittest_$rev.log
	else
		echo "revision still at $rev"
	fi
}

while sleep 1; do
	updaterosetta $HOME/svn/
	updaterosetta $HOME/
	echo OK
	sleep 60
done

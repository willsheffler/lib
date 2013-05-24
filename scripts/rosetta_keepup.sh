#!/bin/bash

# options and globals
CMD_PREF=$1
if [ ! $CMD_PREF ]; then CMD_PREF="nice -n +99"; fi
echo "LOG CMD_PREF: $CMD_PREF"

CMD_SUFF=$2
if [ ! $CMD_SUFF ]; then CMD_SUFF=''; fi
echo "LOG CMD_SUFF: $CMD_SUFF"

SVNURL="https://svn.rosettacommons.org/source/trunk/rosetta"

function ncpu {
	if [ `uname` = "Darwin" ]; then
		ncpu=$(sysctl hw.ncpu | awk '{print $2}')
	else
		ncpu=$(cat /proc/cpuinfo | grep processor | wc -l | grep -v 0)
	fi
	load=$(uptime | egrep -o '\b[0-9]+[.][0-9]+$'|egrep -o '^[0-9]+')
	diff=$[ $ncpu - $load ]
	echo -e "1\n$diff" | sort -nr | head -n1
}

function getsvnrev {
	echo $(svn info $1 | grep 'Last Changed Rev:' | sed -e 's=Last Changed Rev: ==g')
}



function updaterosetta {

	echo "LOG START_UPDATE_ROSETTA_ALL $(date)"

	# get server svn revisions

	revtop=$(getsvnrev $SVNURL                 ); echo "LOG LATEST REV $revtop for $SVNURL"
	revsrc=$(getsvnrev $SVNURL/rosetta_source  ); echo "LOG LATEST REV $revsrc for $SVNURL/rosetta_source"
	revtst=$(getsvnrev $SVNURL/rosetta_tests   ); echo "LOG LATEST REV $revtst for $SVNURL/rosetta_tests"
	revdat=$(getsvnrev $SVNURL/rosetta_database); echo "LOG LATEST REV $revdat for $SVNURL/rosetta_database"

	# loop over src dirs to check for updates

	for ROSETTA_DIR in $@; do

		LBL=$(basename $ROSETTA_DIR)
		L1="LOG     $LBL:"
		R1="RUN     $LBL:"
		L2="LOG         $LBL:"
		R2="RUN         $LBL:"
		L3="LOG             $LBL:"
		R3="RUN             $LBL:"

		echo "$L1 START_UPDATE_ROSETTA_$LBL $(date)"
	

		D="$ROSETTA_DIR/rosetta"
		if [ ! -d $D ]; then echo "ERR ERROR DIRECTORY $D DOES NOT EXIST!"; return -1; fi

		LOG="$D/keepup_logs/${LBL}_r${revtop}_"
		echo "$R1 mkdir -p $D/keepup_logs"
		          mkdir -p $D/keepup_logs

		# get local svn revisions
		localrevtop=1 # $(getsvnrev $D)
		localrevsrc=1 # $(getsvnrev $D/rosetta_source  )
		localrevtst=1 # $(getsvnrev $D/rosetta_tests   )
		localrevdat=1 # $(getsvnrev $D/rosetta_database)

		if [ $revtop = $localrevtop ]; then
			echo "$L2 NO UPDATES NECESSARY: rosetta STILL AT REVISION r$revtop"
		else

			echo "$L2 UPDATING SVN: rosetta UPDATED FROM r$localrevtop TO r$revtop"
			echo "$R2 svn up -r $revtop $D > ${LOG}0_svn_up.log 2> ${LOG}0_svn_up.err"
			$CMD_PREF svn up -r $revtop $D > ${LOG}0_svn_up.log 2> ${LOG}0_svn_up.err $CMD_SUFF

			# now compile if necessary (could be rosetta_source not updated)
			echo "$L2 CHECKING BUILD STATUS"

			if [ $revsrc = $localrevsrc ]; then
				echo "$L3 NOT BUILDING: rosetta_source STILL AT REVISION r$revsrc"
			else
				echo "$L3 REBUILDING AND UNITTESTING: rosetta_source UPDATED FROM r$localrevsrc TO r$revsrc"

				echo "$R3 cd $D/rosetta_source"
				          cd $D/rosetta_source

				echo "$L3 NOT UPDATING OPTIONS, KEEP CURRENT OPTION .gen.hh FILES"
				# echo "$R3 ./update_options.sh > ${LOG}1_update_options.log 2> ${LOG}1_update_options.err"
				# $CMD_PREF ./update_options.sh > ${LOG}1_update_options.log 2> ${LOG}1_update_options.err $CMD_SUFF

				echo "$R3 ./scons.py -kj$(ncpu) mode=release bin            > ${LOG}2_build_release.log     2> ${LOG}2_build_release.err"
				$CMD_PREF ./scons.py -kj$(ncpu) mode=release bin            > ${LOG}2_build_release.log     2> ${LOG}2_build_release.err      $CMD_SUFF
				echo "$R3 ./scons.py -kj$(ncpu) mode=release extras=omp bin > ${LOG}3_build_release_omp.log 2> ${LOG}3_build_release_omp.err"
				$CMD_PREF ./scons.py -kj$(ncpu) mode=release extras=omp bin > ${LOG}3_build_release_omp.log 2> ${LOG}3_build_release_omp.err  $CMD_SUFF
				echo "$R3 ./scons.py -kj$(ncpu) mode=debug                  > ${LOG}4_build_debug.log       2> ${LOG}4_build_debug.err"
				$CMD_PREF ./scons.py -kj$(ncpu) mode=debug                  > ${LOG}4_build_debug.log       2> ${LOG}4_build_debug.err        $CMD_SUFF
				echo "$R3 ./scons.py -kj$(ncpu) mode=debug cxx=clang bin    > ${LOG}5_build_debug_clang.log 2> ${LOG}5_build_debug_clang.err"
				$CMD_PREF ./scons.py -kj$(ncpu) mode=debug cxx=clang bin    > ${LOG}5_build_debug_clang.log 2> ${LOG}5_build_debug_clang.err  $CMD_SUFF
				echo "$R3 ./scons.py -kj$(ncpu) cat=test                    > ${LOG}6_build_test.log        2> ${LOG}6_build_test.err"
				$CMD_PREF ./scons.py -kj$(ncpu) cat=test                    > ${LOG}6_build_test.log        2> ${LOG}6_build_test.err         $CMD_SUFF

				RELEASEOK=$(tail -n1 ${LOG}2_build_release.log)
				DEBUGOK=$(tail -n1 ${LOG}4_build_debug.log)
				TESTOK=$(tail -n1 ${LOG}6_build_test.log)								
				
				echo "RELEASEOK $RELEASEOK, DEBUGOK $DEBUGOK, TESTOK $TESTOK"

				echo "$R3 python test/run.py --database $HOME/tmp/tmp_rosetta/rosetta_database -j $(ncpu) > ${LOG}6_test_unittest.log 2> ${LOG}6_test_unittest.err"
				$CMD_PREF python test/run.py --database $HOME/tmp/tmp_rosetta/rosetta_database -j $(ncpu) > ${LOG}6_test_unittest.log 2> ${LOG}6_test_unittest.err $CMD_SUFF

			fi # if [ $revsrc = $localrevsrc ]

			# now run integration tests if necessary (could be update to tools or demos)
			echo "$L2 CHECKING INTEGRATION TEST STATUS"

			if [ $revsrc = $localrevsrc -a $revdat = $localrevdat -a $revtst = $localrevtst ]; then
				echo "$L3 NOT RUNNING INTEGRATION TESTS: source, database, AND tests STILL AT REVISIONS r$localrevsrc, r$localrevdat, AND r$localrevtst"
			else

				if [ $revsrc != $localrevsrc ]; then echo "$L3 RUNNING INTEGRATION TESTS: rosetta_source UPDATED FROM r$localrevsrc TO r$revsrc"; fi
				if [ $revdat != $localrevdat ]; then echo "$L3 RUNNING INTEGRATION TESTS: rosetta_database UPDATED FROM r$localrevdat TO r$revdat"; fi
				if [ $revtst != $localrevtst ]; then echo "$L3 RUNNING INTEGRATION TESTS: rosetta_tests UPDATED FROM r$localrevtst TO r$revtst"; fi

				echo "$L3 FOR CONSISTENCY IN TEST RESULTS, LINK TO $HOME/tmp/tmp_rosetta FOR TESTING"
				echo "$R3 mkdir -p $HOME/tmp; touch $HOME/tmp/tmp_rosetta; rm -f $HOME/tmp/tmp_rosetta; ln -s $D $HOME/tmp/tmp_rosetta"
				          mkdir -p $HOME/tmp; touch $HOME/tmp/tmp_rosetta; rm -f $HOME/tmp/tmp_rosetta; ln -s $D $HOME/tmp/tmp_rosetta

				echo "$R3 cd $HOME/tmp/tmp_rosetta/rosetta_tests/integration"
				          cd $HOME/tmp/tmp_rosetta/rosetta_tests/integration

				echo "$R3 python integration.py --database $HOME/tmp/tmp_rosetta/rosetta_database -j $(ncpu) > ${LOG}7_integration_test.log 2> ${LOG}7_integration_test.err"
				$CMD_PREF python integration.py --database $HOME/tmp/tmp_rosetta/rosetta_database -j $(ncpu) > ${LOG}7_integration_test.log 2> ${LOG}7_integration_test.err $CMD_SUFF

				echo "$L3 ARCHIVE RESULTS TO: $D/rosetta_tests/integration/integration_test_results_${LBL}_r$revtop.tar.bz2"
				echo "$R3 tar cjf integration_test_results_${LBL}_r$revtop.tar.bz2 new > ${LOG}8_archive_test_results.log 2> ${LOG}8_archive_test_results.err"
				$CMD_PREF tar cjf integration_test_results_${LBL}_r$revtop.tar.bz2 new > ${LOG}8_archive_test_results.log 2> ${LOG}8_archive_test_results.err $CMD_SUFF

			fi # if [ $revsrc = $localrevsrc -a $revdat = $localrevdat -a $revtst = $localrevtst ]

			echo "$R2 find $(dirname $LOG) -empty -type f -print0 | xargs -0 rm"
			nice      find $(dirname $LOG) -empty -type f -print0 | xargs -0 rm

			echo "$L1 DONE_UPDATE_ROSETTA_$LBL $(date)"

		fi # if [ $revtop = $localrevtop ]

	done # for ROSETTA_DIR in $@

	echo "LOG DONE_UPDATE_ROSETTA_ALL $(date)"
}

updaterosetta $HOME/dev $HOME/svn

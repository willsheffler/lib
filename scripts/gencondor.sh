function gencondor {
    for i in 1b72 1di2 1dtj 1hz6 1mky 1n0u 1ogw 2tif; do 
	sed -e "s/\*\*\*\*/$i/g" $1 > condor/$1.$i 
    done
}

function gendfrun {
    for i in `/bin/ls $2`; do 
	file=`basename $i`
	echo "sed -e "s/OUTFILE/$file/g" $1 > condor/$1.$file"
	sed -e "s/OUTFILE/$file/g" $1 > condor/$1.$file 
    done    
}

function genstats {
    template=$1
    date=$2
    lists=$3
    args=$4
    for list in `/bin/ls $lists`; do 
	file=${date}-genstats-`basename $template`-`basename $list`
	fasc=$file
	#echo $date
	#echo $list
	#echo $args
	#echo $file
	sed -e "s:LIST:$list:g" -e "s:LIST:$list:g" -e "s:FASC:$file:g" -e "s:LOGFILE:$file:g" -e "s:ARGS:$args:g" $template > condor/$file.condor
    done
}


function gencsa {
    for i in `/bin/ls prots`; 
    
}

function genpackstats {
	template=$1
	date=$2
	list=$3
	prot=$4
	chain=$5
	sed -e "s:LIST:$list:g" -e "s:DATE:$date:g" -e "s:PROT:$prot:g" -e "s:CHAIN:$chain:g" $template 
	

}


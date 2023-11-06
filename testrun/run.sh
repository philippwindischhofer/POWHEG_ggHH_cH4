#!/bin/bash

# example of a run script for a parallel run on ncores cpu's
ncores=1

# number of iterations for the calculation of xgrid at parallelstage 1. It is the old ncall1
nxgriditeration=3

if [ "$1" = "help" ] || [ "$1" = "" ]
then
    echo "*************************************************"
    echo "-- Menu for running powheg --                    "
    echo "Possible running options:                        "
    echo "./run 0 : run in Higgs Effective Theory (HEFT)   "
    echo "./run 1 : run in Born improved HEFT              "
    echo "./run 2 : run in approx. full theory (FTapp)     "
    echo "./run 3 : run in full theory (grids for virtual) "
    echo "*************************************************"
    exit
else
    mode=$1
fi

function runthem {
    export PYTHONPATH=$PWD:$PYTHONPATH
    for i in $(seq 1 $ncores)
    do
	if [ "$1" = "3" ]
	then
	    if [ ! -f events.cdf ]
	    then
		ln -s ../Virtual/events.cdf events.cdf
	    fi
	    if [ ! -f creategrid.py ]
	    then
		ln -s ../Virtual/creategrid.py creategrid.py
	    fi
	    if [ "`echo Virt_full_*E*.grid`" == "Virt_full_*E*.grid" ]
	    then
		for grid in ../Virtual/Virt_full_*E*.grid; do ln -s $grid ${grid##*/}; done
	    fi
	fi
	# launch powheg
	echo $i | $prg > $logfile-$i.log 2>&1 &
    done
}

prg=../pwhg_main

begin=$(date +"%s")

parstage=1
echo "***********************************************"
echo " stage " $parstage
echo "***********************************************"

for xgrid in $(seq 1 $nxgriditeration)
do
    cat powheg.input-save | sed "s/mtdep.*/mtdep $mode/ ; s/parallelstage.*/parallelstage $parstage/ ; s/xgriditeration.*/xgriditeration $xgrid/">powheg.input
    cp powheg.input powheg.input-$parstage-$xgrid
    logfile=run-$parstage-$xgrid
    runthem $mode
    wait
    rm powheg.input
done

parstage=2
echo "***********************************************"
echo " stage " $parstage
echo "***********************************************"

xgrid=$nxgriditeration

cat powheg.input-save | sed "s/mtdep.*/mtdep $mode/ ; s/parallelstage.*/parallelstage $parstage/ ; s/xgriditeration.*/xgriditeration $xgrid/">powheg.input
cp powheg.input powheg.input-$parstage
logfile=run-$parstage
runthem $mode
wait
rm powheg.input

# termin=$(date +"%s")
# difftimelps=$(($termin-$begin))
# echo "$(($difftimelps / 60)) minutes and $(($difftimelps % 60)) seconds elapsed since script start."
# #exit

# parstage=3
# echo "***********************************************"
# echo " stage " $parstage
# echo "***********************************************"

# cat powheg.input-save | sed "s/mtdep.*/mtdep $mode/ ; s/parallelstage.*/parallelstage $parstage/">powheg.input
# cp powheg.input powheg.input-$parstage
# logfile=run-$parstage
# runthem $mode
# wait
# rm powheg.input


# parstage=4
# echo "***********************************************"
# echo " stage " $parstage
# echo "***********************************************"

# cat powheg.input-save | sed "s/mtdep.*/mtdep $mode/ ; s/parallelstage.*/parallelstage $parstage/ ; s/testplots.*/#testplots 1/">powheg.input
# cp powheg.input powheg.input-$parstage
# logfile=run-$parstage
# runthem $mode
# wait

# termin=$(date +"%s")
# difftimelps=$(($termin-$begin))
# echo "$(($difftimelps / 60)) minutes and $(($difftimelps % 60)) seconds elapsed since script start."

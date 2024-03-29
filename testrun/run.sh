#!/bin/bash

# example of a run script for a parallel run on ncores cpu's
ncores=8

# number of iterations for the calculation of xgrid at parallelstage 1. It is the old ncall1
nxgriditeration=3

if [ "$1" = "help" ]
then
    echo "******************************************************"
    echo "-- Menu for running powheg --                         "
    echo "Possible running options:                             "
    echo "./run warmup : just run warmup phase (generates links/grid)"
    echo "./run all    : run everything                         "
    echo "*************************************************"
    exit
else
    mode=$1
fi

function warmup {
    export PYTHONPATH=$PWD:$PYTHONPATH
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

    # get formatted coupling values
    chhh=$(awk 'sub(/^chhh/,""){printf "%+.4E", $1}' powheg.input-save)
    ct=$(printf "%+.4E" 1.0)
    ctt=$(printf "%+.4E" 0.0)
    cg=$(printf "%+.4E" 0.0)
    cgg=$(printf "%+.4E" 0.0)

    gridtemp="Virt_full_${chhh}_${ct}_${ctt}_${cg}_${cgg}.grid"

    pythoncmd="import creategrid as cg; cg.combinegrids('${gridtemp}', ${chhh}, ${ct}, ${ctt}, ${cg}, ${cgg})"
    python3 -c "$pythoncmd"
}

function runthem {
    for i in $(seq 1 $ncores)
    do
	# launch powheg
	echo $i | $prg > $logfile-$i.log 2>&1 &
    done
}

prg=../pwhg_main

begin=$(date +"%s")

if [ "$mode" = "warmup" ] || [ "$mode" = "all" ]
then
    echo "***********************************************"
    echo " stage warmup"
    echo "***********************************************"

    warmup
    wait

    if [ "$mode" = "warmup" ]
    then
        exit
    fi
fi

parstage=1
echo "***********************************************"
echo " stage " $parstage
echo "***********************************************"

for xgrid in $(seq 1 $nxgriditeration)
do
    cat powheg.input-save | sed "s/parallelstage.*/parallelstage $parstage/ ; s/xgriditeration.*/xgriditeration $xgrid/">powheg.input
    cp powheg.input powheg.input-$parstage-$xgrid
    logfile=run-$parstage-$xgrid
    runthem
    wait
    rm powheg.input
done


parstage=2
echo "***********************************************"
echo " stage " $parstage
echo "***********************************************"

cat powheg.input-save | sed "s/parallelstage.*/parallelstage $parstage/ ; s/xgriditeration.*/xgriditeration $xgrid/">powheg.input
cp powheg.input powheg.input-$parstage
logfile=run-$parstage
runthem
wait
rm powheg.input

termin=$(date +"%s")
difftimelps=$(($termin-$begin))
echo "$(($difftimelps / 60)) minutes and $(($difftimelps % 60)) seconds elapsed since script start."
#exit

parstage=3
echo "***********************************************"
echo " stage " $parstage
echo "***********************************************"

cat powheg.input-save | sed "s/parallelstage.*/parallelstage $parstage/ ; s/xgriditeration.*/xgriditeration $xgrid/">powheg.input
cp powheg.input powheg.input-$parstage
logfile=run-$parstage
runthem
wait
rm powheg.input


parstage=4
echo "***********************************************"
echo " stage " $parstage
echo "***********************************************"

cat powheg.input-save | sed "s/parallelstage.*/parallelstage $parstage/ ; s/xgriditeration.*/xgriditeration $xgrid/ ; s/testplots.*/#testplots 1/">powheg.input
cp powheg.input powheg.input-$parstage
logfile=run-$parstage
runthem
wait

termin=$(date +"%s")
difftimelps=$(($termin-$begin))
echo "$(($difftimelps / 60)) minutes and $(($difftimelps % 60)) seconds elapsed since script start."

#!/bin/bash

# example of a run script for a parallel run on ncores cpu's
ncores=8

# number of iterations for the calculation of xgrid at parallelstage 1. It is the old ncall1
nxgriditeration=3

if [ "$1" = "help" ] || [ "$1" = "" ]
then
    echo "******************************************************"
    echo "-- Menu for running powheg --                         "
    echo "Possible running options:                             "
    echo "./run warmup : run warmup phase (generates links/grid)"
    echo "./run 0      : run in Higgs Effective Theory (HEFT)   "
    echo "./run 1      : run in Born improved HEFT              "
    echo "./run 2      : run in approx. full theory (FTapp)     "
    echo "./run 3      : run in full theory (grids for virtual) "
    echo "*************************************************"
    exit
else
    mode=$1
fi

function warmup {
    export PYTHONPATH=$PWD:$PYTHONPATH
    # creates symbolic links to the files that are needed
    # for mtdep = 3
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
    ct=$(awk 'sub(/^ct /,""){printf "%+.4E", $1}' powheg.input-save)
    ctt=$(awk 'sub(/^ctt/,""){printf "%+.4E", $1}' powheg.input-save)
    cg=$(awk 'sub(/^cggh /,""){printf "%+.4E", $1}' powheg.input-save)
    cgg=$(awk 'sub(/^cgghh/,""){printf "%+.4E", $1}' powheg.input-save)

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

if [ "$mode" = "warmup" ] || [ "$mode" = "3" ]
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
    cat powheg.input-save | sed "s/mtdep.*/mtdep $mode/ ; s/parallelstage.*/parallelstage $parstage/ ; s/xgriditeration.*/xgriditeration $xgrid/">powheg.input
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

cat powheg.input-save | sed "s/mtdep.*/mtdep $mode/ ; s/parallelstage.*/parallelstage $parstage/ ; s/xgriditeration.*/xgriditeration $xgrid/">powheg.input
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

cat powheg.input-save | sed "s/mtdep.*/mtdep $mode/ ; s/parallelstage.*/parallelstage $parstage/ ; s/xgriditeration.*/xgriditeration $xgrid/">powheg.input
cp powheg.input powheg.input-$parstage
logfile=run-$parstage
runthem
wait
rm powheg.input


parstage=4
echo "***********************************************"
echo " stage " $parstage
echo "***********************************************"

cat powheg.input-save | sed "s/mtdep.*/mtdep $mode/ ; s/parallelstage.*/parallelstage $parstage/ ; s/xgriditeration.*/xgriditeration $xgrid/ ; s/testplots.*/#testplots 1/">powheg.input
cp powheg.input powheg.input-$parstage
logfile=run-$parstage
runthem
wait

termin=$(date +"%s")
difftimelps=$(($termin-$begin))
echo "$(($difftimelps / 60)) minutes and $(($difftimelps % 60)) seconds elapsed since script start."

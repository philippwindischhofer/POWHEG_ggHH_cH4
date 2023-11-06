#!/bin/bash

# Store commandline arguments
STARTSEED="$1"
PRG="$2"
LOGFILE="$3"
echo "[${SLURM_PROCID} (`hostname`)] Received: STARTSEED=${STARTSEED}, PRG=${PRG}, LOGFILE=${LOGFILE}"

# Compute seed
SEED=$((${STARTSEED} + ${SLURM_PROCID}))

# Run pwhg
echo "[${SLURM_PROCID} (`hostname`)] Running: echo ${SEED} | ${PRG} > ${LOGFILE}-${SEED}.log 2>&1"
echo ${SEED} | ${PRG} > ${LOGFILE}-${SEED}.log 2>&1

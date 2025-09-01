#!/bin/bash
# Define common options

# === Detect cluster and set account/partition ===
ENV=""
if [[ "$PWD" == "/leonardo/home/userexternal/gbillo00"* ]]; then
	COMMON_OPTS="--time=02:00:00 \
	     	     --exclusive
                     --partition=dcgp_usr_prod 
		     -A uTS25_Tornator_0"
  ENV="LEONARDO"
  module load gcc/12.2.0
  module load openmpi/4.1.6--gcc--12.2.0
  echo "Detected Leonardo environment"
elif [[ "$PWD" == "/u/dssc/gbillo/HPCproject"* ]]; then
  COMMON_OPTS=" --time=02:00:00 \
	     	--exclusive
  		--partition=EPYC 
		-A dssc"
  ENV="ORFEO"
  echo "Detected Orfeo environment"
  module load openMPI
else
  echo "❌ Unknown system environment (PWD=$PWD)" >&2
  exit 1
fi
# COMMON_OPTS=" --time=02:00:00 \
# 	     --exclusive"

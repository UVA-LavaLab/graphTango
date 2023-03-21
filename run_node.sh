#!/bin/bash

if [ $# -ne 0 ]; then
	echo "Usage: ./run_slurm.sh <node>"
	exit -1
fi

SLURM_TEMPLATE=slurm_node.template

export FULCRUM_NODE=$1

function run_config {

	OUT_DIR="gapbs_outputs"
	mkdir -p $OUT_DIR
	

	#copy the slurm script
	envsubst < $SLURM_TEMPLATE > $OUT_DIR/slurm.sh

	#add to queue
	sbatch --chdir $OUT_DIR slurm.sh
}

#run_config 64
#run_config 128
#run_config 256
#run_config 512
#run_config 1024
#run_config 2048
#run_config 4096
#run_config 8192
#run_config 16384
#run_config 32768
run_config


#!/bin/bash
# As precaution, clear OMP stuff and remove any alg or update csv files in the folder
unset OMP_DISPLAY_ENV OMP_NUM_THREADS OMP_PROC_BIND OMP_PLACES
rm -rf Alg*.csv
rm -rf Update*.csv 

#export LD_LIBRARY_PATH=/home/alif/.installs/likwid/lib/ib/ 
export OMP_DISPLAY_ENV=true
export OMP_NUM_THREADS=12
export OMP_PROC_BIND=close
#export OMP_PLACES={2}:64:1
export OMP_PLACES=threads

dataDir=../../
#STRUCTURES=(adListChunked adListShared degAwareRHH stinger)
#STRUCTURES=(graphite adListChunked adListShared degAwareRHH stinger)
STRUCTURES=(graphTango)
batchSize=1000000
NumThreads=12

# whether each algorithm is weighted or unweighted
declare -A ALGORITHMS 
ALGORITHMS=(         
#         [ccdyn]=0         
#         [prdyn]=0
         [bfsdyn]=0
#         [mcdyn]=0 
#         [sswpdyn]=1
 #        [ssspdyn]=1
)

# Max num_nodes to initialize for each dataset
declare -A DATASETS
DATASETS=(
       [orkut.el]=3072441
)

  for dataset in "${!DATASETS[@]}"; do  
    for structure in "${STRUCTURES[@]}"; do
      for algorithm in "${!ALGORITHMS[@]}"; do 
         echo ../frontEnd -d 1 -w ${ALGORITHMS[$algorithm]} -f ${dataDir}$dataset -b ${batchSize} -s $structure -n ${DATASETS[$dataset]} -a $algorithm -t ${NumThreads}   
         ./frontEnd -d 1 -o "update.csv" -w ${ALGORITHMS[$algorithm]} -f ${dataDir}$dataset -b ${batchSize} -s $structure -n ${DATASETS[$dataset]} -a $algorithm -t ${NumThreads}  
        DIRECTORY=./results/${algorithm}/${structure}/${dataset}/
        mkdir -p ${DIRECTORY}
        mv Alg*.csv ${DIRECTORY}/
        mv update*.csv ${DIRECTORY}/
      done    
    done  
  done 


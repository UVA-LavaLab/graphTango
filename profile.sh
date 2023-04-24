#!/bin/bash
# As precaution, clear OMP stuff and remove any alg or update csv files in the folder
unset OMP_DISPLAY_ENV OMP_NUM_THREADS OMP_PROC_BIND OMP_PLACES
rm -rf time.csv

#export LD_LIBRARY_PATH=/home/alif/.installs/likwid/lib/ib/ 
export OMP_DISPLAY_ENV=true
export OMP_NUM_THREADS=24
export OMP_PROC_BIND=close
#export OMP_PLACES={2}:64:1
export OMP_PLACES=threads

dataDir=/bigtemp/fas9nw/datasets/
#STRUCTURES=(adListChunked adListShared stinger)
#STRUCTURES=(graphite adListChunked adListShared degAwareRHH stinger)
STRUCTURES=(graphTango)
batchSize=1000000
NumThreads=24

# whether each algorithm is weighted or unweighted
declare -A ALGORITHMS 
ALGORITHMS=(         
         [ccdyn]=0         
         [prdyn]=0
         [bfsdyn]=0
#         [mcdyn]=0 
#         [sswpdyn]=1
         [ssspdyn]=1
)

# Max num_nodes to initialize for each dataset
declare -A DATASETS
DATASETS=(
       [liveJournal.el]=4847571       
       [orkut.el]=3072441       
       [wiki.el]=1791489
       [talk.el]=2394385
#       [twitter.el]=61578415       
       [road.el]=23947347       
#       [web.el]=50636151
       [urand_24.el]=16777216
       [kron_24.el]=16777216
)

  for dataset in "${!DATASETS[@]}"; do  
    for structure in "${STRUCTURES[@]}"; do
      for algorithm in "${!ALGORITHMS[@]}"; do 
         echo ../frontEnd -d 1 -w ${ALGORITHMS[$algorithm]} -f ${dataDir}$dataset -b ${batchSize} -s $structure -n ${DATASETS[$dataset]} -a $algorithm -t ${NumThreads}   
         ../frontEnd -d 1 -w ${ALGORITHMS[$algorithm]} -f ${dataDir}$dataset -b ${batchSize} -s $structure -n ${DATASETS[$dataset]} -a $algorithm -t ${NumThreads}  
        DIRECTORY=./results/${algorithm}/${structure}/${dataset}/
        mkdir -p ${DIRECTORY}
        mv time.csv ${DIRECTORY}/
      done    
    done  
  done 


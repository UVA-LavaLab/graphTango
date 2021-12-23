#!/bin/bash
# As precaution, clear OMP stuff and remove any alg or update csv files in the folder
unset OMP_DISPLAY_ENV OMP_NUM_THREADS OMP_PROC_BIND OMP_PLACES
rm Alg*.csv
rm Update*.csv 

export OMP_DISPLAY_ENV=true
export OMP_NUM_THREADS=16
export OMP_PROC_BIND=close
#export OMP_PLACES={2}:64:1
export OMP_PLACES=threads

dataDir=./inputResource/
#STRUCTURES=(adListChunked adListShared degAwareRHH stinger)
#STRUCTURES=(graphite)
#STRUCTURES=(stinger)
STRUCTURES=(stinger)
batchSize=1000000 
RUNS=3
NumThreads=16

# whether each algorithm is weighted or unweighted
declare -A ALGORITHMS 
ALGORITHMS=(         
#         [ccdyn]=0         
#         [ccfromscratch]=0 
#         [prdyn]=0
#         [prfromscratch]=0         
         [bfsdyn]=0
#         [bfsfromscratch]=0 
#         [mcdyn]=0 
#         [mcfromscratch]=0
#         [sswpdyn]=1
#         [sswpfromscratch]=1
#         [ssspdyn]=1
#         [ssspfromscratch]=1         
)

# Max num_nodes to initialize for each dataset
declare -A DATASETS
DATASETS=(
       [LiveJournal.csv]=4847571       
       [orkut.csv]=3072441       
       [wiki-topcats.csv]=1791489
       [WikiTalk.csv]=2394385
)

runs=${RUNS}
while [ $runs -gt 0 ]
do 
  for dataset in "${!DATASETS[@]}"; do  
    for structure in "${STRUCTURES[@]}"; do
      for algorithm in "${!ALGORITHMS[@]}"; do 
        # if-else to make sure orkut runs in undirected mode 
        #if [ "$dataset" == "com-orkut.ungraph.shuffle.t.w.csv" ]; 
        #then
        # echo ./frontEnd -d 0 -w ${ALGORITHMS[$algorithm]} -f ${dataDir}$dataset -b ${batchSize} -s $structure -n ${DATASETS[$dataset]} -a $algorithm -t ${NumThreads}   
        # ./frontEnd -d 0 -w ${ALGORITHMS[$algorithm]} -f ${dataDir}$dataset -b ${batchSize} -s $structure -n ${DATASETS[$dataset]} -a $algorithm -t ${NumThreads}   
        #else 
         echo ./frontEnd -d 1 -w ${ALGORITHMS[$algorithm]} -f ${dataDir}$dataset -b ${batchSize} -s $structure -n ${DATASETS[$dataset]} -a $algorithm -t ${NumThreads}   
         ./frontEnd -d 1 -w ${ALGORITHMS[$algorithm]} -f ${dataDir}$dataset -b ${batchSize} -s $structure -n ${DATASETS[$dataset]} -a $algorithm -t ${NumThreads}
        #fi        
        # make right directory and move the two generated files into it  
        DIRECTORY=./results/${algorithm}/${structure}/${dataset}/Run${runs}      
        mkdir -p ${DIRECTORY}
        mv Alg*.csv ${DIRECTORY}/
        mv Update*.csv ${DIRECTORY}/
      done    
    done  
  done 
runs=$(( $runs - 1 ))
done 

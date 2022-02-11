#!/bin/bash
# As precaution, clear OMP stuff and remove any alg or update csv files in the folder
unset OMP_DISPLAY_ENV OMP_NUM_THREADS OMP_PROC_BIND OMP_PLACES
rm Alg*.csv
rm Update*.csv 

export LD_LIBRARY_PATH=/home/alif/.installs/likwid/lib/ib/ 
export OMP_DISPLAY_ENV=true
export OMP_NUM_THREADS=12
export OMP_PROC_BIND=close
#export OMP_PLACES={2}:64:1
export OMP_PLACES=threads

dataDir=./inputResource/
#STRUCTURES=(adListChunked adListShared degAwareRHH stinger)
#STRUCTURES=(graphite adListChunked adListShared degAwareRHH stinger)
STRUCTURES=(degAwareRHH)
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
#         [ssspdyn]=1
)

# Max num_nodes to initialize for each dataset
declare -A DATASETS
DATASETS=(
       [LiveJournal.csv]=4847571       
       [orkut.csv]=3072441       
       [wiki-topcats.csv]=1791489
       [WikiTalk.csv]=2394385
#		[rmat_1_1.csv]=1048576
#		[rmat_1_2.csv]=1048576
#		[rmat_1_4.csv]=1048576
#		[rmat_1_8.csv]=1048576
#		[rmat_1_16.csv]=1048576
#		[rmat_1_32.csv]=1048576
#		[rmat_1_64.csv]=1048576
#		[rmat_1_128.csv]=1048576
#		[rmat_1_256.csv]=1048576
#       [rmat.csv]=33554432
)

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
        DIRECTORY=./results/${algorithm}/${structure}/${dataset}/
        mkdir -p ${DIRECTORY}
        mv Alg*.csv ${DIRECTORY}/
        mv Update*.csv ${DIRECTORY}/
      done    
    done  
  done 


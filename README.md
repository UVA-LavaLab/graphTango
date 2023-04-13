
This is the SAGA-Bench integration of **GraphTango**. GraphTango is a hybrid representation format that provides excellent update and analytics throughput regardless of the graph's degree distribution. GraphTango dynamically switches among three different formats based on a vertex's degree: i) Low-degree vertices store the edges directly with the neighborhood metadata, confining accesses to a single cache line, ii)  Medium-degree vertices use adjacency lists, and iii) High-degree vertices use hash tables as well as adjacency lists. In this case, the adjacency list provides fast traversal during the analytics phase, while the hash table provides constant-time lookups during the update phase. Performance is further optimized by designing an open-addressing-based hash table that fully utilizes every fetched cache line.


## Overview of the Directory Structure
Please refer to [SAGA-Bench](https://github.com/abasak24/SAGA-Bench) for the common directory structure. **GraphTango** specific files are the following:
1. **src/dynamic/GraphTango.h**: Contains the implementation of the GraphTango API. Supports insertion/deletion of edges and vertices, both individual and batched.
2. **src/dynamic/Vertex.h**: Implements the data structure partaining to a single vertex.
3. **src/dynamic/LockFreePoolWithList.h**: Custom memory allocator optimized for GraphTango.
4. **src/dynamic/common.h** Contains various configurations of GraphTango (e.g., different hashing mechanism or memory allocators). Default is the cache-friendly-hashing scheme and a custom memory allocator.


## Input Datasets
We used *.csv* format where each line contains the following:
```
[source vertex ID], [destination vertex ID], [timestamp], [weight]
```
Graph datasets are first randomly shuffled to break any ordering in the input files. This is done to ensure the realistic scenario that streaming edges are not likely to come in any pre-defined order. The shuffled input file is then read in batches. The resources for preparing the input datasets are provided in the folder inputResource. `inputResource/shuffle.sh` can be used to shuffle a dataset file in .txt format (e.g., those found in [SNAP](https://snap.stanford.edu/data/)). After shuffling, timestamps and weights can be added using `inputResource/addWeightAndTime.sh` and `inputResource/appendValues.py`, which will result in the final *.csv* format.


## Compiling and Running GraphTango
GraphTango has been tested on Ubuntu 20.04 LTS with gcc 9.3.0. To build, run the following commands:

```
$ git clone https://github.com/alifahmed/graphTango.git
$ cd graphTango/
$ make 
```

An executable `frontEnd` will be created. `frontEnd` should be run with the following parameters. `./frontEnd --help` also provides this information.

```
-f : provides a location to an input graph file in .csv format
-b : batch size
-d : whether the input graph is directed or undirected. 0=undirected; 1=directed.
-w : whether weights should be read from the input file. 0=don't read weights; 1=read weights. Weights are required only for SSSP and SSWP. 
-s : data structure to be used (see DATA STRUCTURE OPTIONS below). 
-a : algorithm to be run (see ALGORITHM OPTIONS below). 
-n : max number of nodes the data structure must be initialized with. 
-t : max number of allowed threads.

DATA STRUCTURE OPTIONS: 1) adListShared 2) adListChunked 3) degAwareRHH 4) stinger 5) graphTango (default)
ALGORITHM OPTIONS: 1) prfromscratch 2) prdyn 3) ccfromscratch 4) ccdyn 5) mcfromscratch 6) mcdyn 7) bfsfromscratch 8) bfsdyn (default) 9) ssspfromscratch 10) ssspdyn 11) sswpfromscratch 12) sswpdyn
```

For example, to run BFS using GraphTango, the following command can be used:
```
$ ./frontEnd -f example.csv -b 1000000 -w 0 -d 1 -s graphTango -n 5000000 -a bfsdyn -t 16
```

Each run generates two csv files: **Alg.csv** and **Update.csv**. These files contain per-batch analytics and update times, respectively, in seconds.

## How to Cite

If you are using this, please cite:
Alif Ahmed, Farzana Ahmed Siddique, Kevin Skadron. [*GraphTango: A Hybrid Representation Format for Efficient Streaming Graph Updates and Analysis*](https://arxiv.org/abs/2212.11935). arXiv:2212.11935 [cs.DS], 2022.


## Contact
In case of any concern, please contact Alif Ahmed at alifahmed@virginia.edu.

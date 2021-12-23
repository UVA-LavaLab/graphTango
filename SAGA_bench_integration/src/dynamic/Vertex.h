#pragma once

#include <cstdlib>
#include <vector>
#include "common.h"
#include "CustomAllocator.h"
#include "omp.h"

typedef std::unordered_map<Idx, u64, std::unordered_map<Idx, u64>::hasher, std::unordered_map<Idx, u64>::key_equal, custom_allocator< std::pair<const Idx,u64>> > graphite_hashmap;
//typedef std::unordered_map<Idx, u64> graphite_hashmap;


#ifdef USE_HYBRID_HASHMAP_WITH_GROUPING

template <typename Neigh>
class EdgeArray{
public:
	u64 							degree = 0;
	u64 							capacity = 0;
	Neigh* 				__restrict 	neighArr = nullptr;
	graphite_hashmap* 	__restrict 	mapArr = nullptr;
};


template <typename Neigh>
class Vertex{
public:
	 EdgeArray<Neigh>		inEdges;
	 EdgeArray<Neigh>		outEdges;
};


#endif


#ifdef USE_HYBRID_HASHMAP_WITH_GROUPING_TIGHTER

#define INITIAL_EDGES				5

template <typename Neigh>
class EdgeArray{
public:
	u32 							degree = 0;
	u32 							capacity = INITIAL_EDGES;
	Neigh* 				__restrict 	neighArr = nullptr;
	graphite_hashmap* 	__restrict 	mapArr = nullptr;
	Neigh							neigh[INITIAL_EDGES];
};


template <typename Neigh>
class Vertex{
public:
	 EdgeArray<Neigh>		inEdges;
	 EdgeArray<Neigh>		outEdges;
};


#endif


#ifdef USE_HYBRID_HASHMAP_WITH_GROUPING_AND_EDGE_ARR_LOCKING

template <typename Neigh>
class EdgeArray{
public:
	u32 							degree = 0;
	u32 							capacity = 0;
	Neigh* 				__restrict 	neighArr = nullptr;
	graphite_hashmap* 	__restrict 	mapArr = nullptr;
	omp_lock_t						lock;
};


template <typename Neigh>
class Vertex{
public:
	 EdgeArray<Neigh>		inEdges;
	 EdgeArray<Neigh>		outEdges;
};

#endif


#ifdef USE_SORTED_EDGES

template <typename Neigh>
class EdgeArray{
public:
	u32 							degree = 0;
	u32 							capacity = 0;
	Neigh* 				__restrict 	neighArr = nullptr;
	graphite_hashmap* 	__restrict 	mapArr = nullptr;
};


template <typename Neigh>
class Vertex{
public:
	 EdgeArray<Neigh>		inEdges;
	 EdgeArray<Neigh>		outEdges;
};

#endif


#ifdef USE_CAHCE_FRIENDLY_HASH

#include "GraphiteHash.h"

template <typename Neigh>
class Vertex{
public:
	 GraphiteHash<Neigh>		inEdges;
	 GraphiteHash<Neigh>		outEdges;
};

#endif


#ifdef USE_CAHCE_FRIENDLY_HASH_ONLY

#include "GraphiteHash.h"

template <typename Neigh>
class Vertex{
public:
	 GraphiteHash<Neigh>		inEdges;
	 GraphiteHash<Neigh>		outEdges;
};

#endif


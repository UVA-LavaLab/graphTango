#pragma once

#include <cstdlib>
#include <vector>
#include "common.h"
#include "GEmpty.h"
#include "CustomAllocator.h"
#include <omp.h>
#include <unordered_map>
#include "Vertex.h"

template <typename Neigh>
class VertexArray {
public:

#ifdef USE_HYBRID_HASHMAP
	u64* 				__restrict outDegree = nullptr;
	u64*				__restrict outCapacity = nullptr;
	Neigh**				__restrict outNeighArr = nullptr;
	graphite_hashmap** 	__restrict outMapArr = nullptr;


	u64* 				__restrict inDegree = nullptr;
	u64*				__restrict inCapacity = nullptr;
	Neigh**				__restrict inNeighArr = nullptr;
	graphite_hashmap** 	__restrict inMapArr = nullptr;

	//u64 vSize = 0;

	void resize(u64 size) {
		//vSize = size;
		outDegree = (u64*)globalAllocator.allocate(size * sizeof(u64));
		outCapacity = (u64*)globalAllocator.allocate(size * sizeof(u64));
		outNeighArr = (Neigh**)globalAllocator.allocate(size * sizeof(Neigh*));
		outMapArr = (graphite_hashmap**)globalAllocator.allocate(size * sizeof(graphite_hashmap*));


		inDegree = (u64*)globalAllocator.allocate(size * sizeof(u64));
		inCapacity = (u64*)globalAllocator.allocate(size * sizeof(u64));
		inNeighArr = (Neigh**)globalAllocator.allocate(size * sizeof(Neigh*));
		inMapArr = (graphite_hashmap**)globalAllocator.allocate(size * sizeof(graphite_hashmap*));
	}

#endif


#ifdef USE_HYBRID_HASHMAP_WITH_CFH
	u64* 				__restrict outDegree = nullptr;
	u64*				__restrict outCapacity = nullptr;
	Neigh**				__restrict outNeighArr = nullptr;
	DstLocPair** 		__restrict outMapArr = nullptr;


	u64* 				__restrict inDegree = nullptr;
	u64*				__restrict inCapacity = nullptr;
	Neigh**				__restrict inNeighArr = nullptr;
	DstLocPair** 		__restrict inMapArr = nullptr;

	//u64 vSize = 0;

	void resize(u64 size) {
		//vSize = size;
		outDegree = (u64*)globalAllocator.allocate(size * sizeof(u64));
		outCapacity = (u64*)globalAllocator.allocate(size * sizeof(u64));
		outNeighArr = (Neigh**)globalAllocator.allocate(size * sizeof(Neigh*));
		outMapArr = (DstLocPair**)globalAllocator.allocate(size * sizeof(DstLocPair*));


		inDegree = (u64*)globalAllocator.allocate(size * sizeof(u64));
		inCapacity = (u64*)globalAllocator.allocate(size * sizeof(u64));
		inNeighArr = (Neigh**)globalAllocator.allocate(size * sizeof(Neigh*));
		inMapArr = (DstLocPair**)globalAllocator.allocate(size * sizeof(DstLocPair*));
	}

#endif


#ifdef USE_ONLY_LINEAR
	u64* 				__restrict outDegree = nullptr;
	u64*				__restrict outCapacity = nullptr;
	Neigh**				__restrict outNeighArr = nullptr;


	u64* 				__restrict inDegree = nullptr;
	u64*				__restrict inCapacity = nullptr;
	Neigh**				__restrict inNeighArr = nullptr;

	//u64 vSize = 0;

	void resize(u64 size) {
		//vSize = size;
		outDegree = (u64*)globalAllocator.allocate(size * sizeof(u64));
		outCapacity = (u64*)globalAllocator.allocate(size * sizeof(u64));
		outNeighArr = (Neigh**)globalAllocator.allocate(size * sizeof(Neigh*));

		inDegree = (u64*)globalAllocator.allocate(size * sizeof(u64));
		inCapacity = (u64*)globalAllocator.allocate(size * sizeof(u64));
		inNeighArr = (Neigh**)globalAllocator.allocate(size * sizeof(Neigh*));
	}
#endif

#ifdef USE_ONLY_HASHMAP
	u64* 				__restrict outDegree = nullptr;
	u64*				__restrict outCapacity = nullptr;
	Neigh**				__restrict outNeighArr = nullptr;
	graphite_hashmap* 	__restrict outMapArr = nullptr;


	u64* 				__restrict inDegree = nullptr;
	u64*				__restrict inCapacity = nullptr;
	Neigh**				__restrict inNeighArr = nullptr;
	graphite_hashmap* 	__restrict inMapArr = nullptr;

	//u64 vSize = 0;

	void resize(u64 size) {
		//vSize = size;
		outDegree = (u64*)globalAllocator.allocate(size * sizeof(u64));
		outCapacity = (u64*)globalAllocator.allocate(size * sizeof(u64));
		outNeighArr = (Neigh**)globalAllocator.allocate(size * sizeof(Neigh*));
		//outMapArr = (graphite_hashmap*)globalAllocator.allocate(size * sizeof(graphite_hashmap*));
		outMapArr = new graphite_hashmap[size];


		inDegree = (u64*)globalAllocator.allocate(size * sizeof(u64));
		inCapacity = (u64*)globalAllocator.allocate(size * sizeof(u64));
		inNeighArr = (Neigh**)globalAllocator.allocate(size * sizeof(Neigh*));
		//inMapArr = (graphite_hashmap**)globalAllocator.allocate(size * sizeof(graphite_hashmap*));
		inMapArr = new graphite_hashmap[size];
	}

#endif

#ifdef USE_WEIRD_SCHEME
	u64* 				__restrict outDegree = nullptr;
	u64*				__restrict outCapacity = nullptr;
	Neigh**				__restrict outNeighArr = nullptr;
	graphite_hashmap* 	__restrict outMapArr = nullptr;


	u64* 				__restrict inDegree = nullptr;
	u64*				__restrict inCapacity = nullptr;
	Neigh**				__restrict inNeighArr = nullptr;
	graphite_hashmap* 	__restrict inMapArr = nullptr;

	//u64 vSize = 0;

	void resize(u64 size) {
		//vSize = size;
		outDegree = (u64*)globalAllocator.allocate(size * sizeof(u64));
		outCapacity = (u64*)globalAllocator.allocate(size * sizeof(u64));
		outNeighArr = (Neigh**)globalAllocator.allocate(size * sizeof(Neigh*));
		//outMapArr = (graphite_hashmap*)globalAllocator.allocate(size * sizeof(graphite_hashmap*));
		outMapArr = new graphite_hashmap[size];


		inDegree = (u64*)globalAllocator.allocate(size * sizeof(u64));
		inCapacity = (u64*)globalAllocator.allocate(size * sizeof(u64));
		inNeighArr = (Neigh**)globalAllocator.allocate(size * sizeof(Neigh*));
		//inMapArr = (graphite_hashmap**)globalAllocator.allocate(size * sizeof(graphite_hashmap*));
		inMapArr = new graphite_hashmap[size];
	}

#endif

};






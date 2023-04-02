#pragma once

#include <cstdlib>
#include <vector>
#include "common.h"
#include "CustomAllocator.h"
#include "omp.h"

#if defined(USE_GT_BALANCED_MALLOC_STDMAP)
typedef std::unordered_map<u32, u32> graphite_hashmap;
#elif defined(USE_GT_BALANCED_MALLOC)

#elif defined(USE_GT_BALANCED_ABSEIL)
#include "absl/container/flat_hash_map.h"
typedef absl::flat_hash_map<u32, u32> graphite_hashmap;

#elif defined(USE_GT_BALANCED_RHH)
#include "robin_hood.h"
typedef robin_hood::unordered_flat_map<u32, u32> graphite_hashmap;

#elif defined(USE_GT_BALANCED_TSL_RHH)
#include "robin_map.h"
typedef tsl::robin_map<u32, u32> graphite_hashmap;

#else
typedef std::unordered_map<u32, u32, std::unordered_map<u32, u32>::hasher, std::unordered_map<u32, u32>::key_equal, custom_allocator< std::pair<const u32,u32>> > graphite_hashmap;
#endif


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


#if defined(USE_GT_BALANCED) || defined(USE_GT_BALANCED_DYN_PARTITION)

#define 	FLAG_EMPTY_SLOT			0xFFFFFFFFU
#define 	FLAG_TOMB_STONE			0xFFFFFFFEU
#define		CACHE_LINE_SIZE			64

template <typename Neigh>
class alignas(CACHE_LINE_SIZE) EdgeArray{

private:
	void rebuildHashTable(u64 oldCap, u64 newCap){
		if(oldCap > TH1){
			//free old map
			globalAllocator.freePow2(etype.type3.mapArr, oldCap * 2 * sizeof(DstLocPair));
			etype.type3.mapArr = nullptr;
		}

		if(newCap > TH1){
			//allocate new map
			etype.type3.mapArr = (DstLocPair*)globalAllocator.allocate(newCap * 2 * sizeof(DstLocPair));

			DstLocPair* __restrict locMap = etype.type3.mapArr;
			memset(locMap, -1, newCap * 2 * sizeof(DstLocPair));

			const u32 mask = newCap * 2 - 1;

			//add existing nodes to hash
			const Neigh* __restrict nn = etype.type3.neighArr;
			for(u64 i = 0; i < degree; i++){
				const u32 dst = nn[i].node;
				u32 idx = dst & mask;
				while(true){
					if(locMap[idx].dst == FLAG_EMPTY_SLOT){
						//found insertion point
						locMap[idx].dst = dst;
						locMap[idx].loc = i;
						break;
					}
					//move on
					idx++;
					if(idx == (newCap * 2)){
						idx = 0;
					}
				}
			}
		}
	}


public:

	const static u64 TH0 = ((CACHE_LINE_SIZE - sizeof(u32) - sizeof(u32)) / sizeof(Neigh));
	const static u64 TH1 = HYBRID_HASH_PARTITION;

	u32 degree = 0;
	u32 capacity = TH0;

	union {
		struct {
			Neigh neigh[TH0];
		} type1;

		struct {
			Neigh* 			__restrict neighArr = nullptr;
		} type2;

		struct {
			Neigh* 			__restrict neighArr = nullptr;
			DstLocPair* 	__restrict mapArr = nullptr;
		} type3;

	} etype;

	u8 __pad[CACHE_LINE_SIZE - sizeof(u32) - sizeof(u32) - sizeof(etype)];

	void insertEdge(const Idx dstId, const Weight weight, u64& edgeCnt){
		//First, check if needs expanding
		if(__builtin_expect(degree == capacity, false)){
			capacity = getNextPow2(capacity * 2);

			if(degree <= TH0){	//Type 1 => Type 2
				Neigh* __restrict newPtr = (Neigh*)globalAllocator.allocPow2(capacity * sizeof(Neigh));
				memcpy(newPtr, etype.type1.neigh, degree * sizeof(Neigh));
				etype.type2.neighArr = newPtr;
			}
			else if(degree <= TH1 && capacity <= TH1) { // Type 2 => Type 2
				Neigh* __restrict newPtr = (Neigh*)globalAllocator.allocPow2(capacity * sizeof(Neigh));
				memcpy(newPtr, etype.type2.neighArr, degree * sizeof(Neigh));
				globalAllocator.freePow2(etype.type2.neighArr, degree * sizeof(Neigh));
				etype.type2.neighArr = newPtr;
			}
			else if(degree <= TH1 && capacity > TH1){ // Type 2 => Type 3
				Neigh* __restrict newPtr = (Neigh*)globalAllocator.allocPow2(capacity * sizeof(Neigh));
				memcpy(newPtr, etype.type2.neighArr, degree * sizeof(Neigh));
				globalAllocator.freePow2(etype.type2.neighArr, degree * sizeof(Neigh));
				etype.type3.neighArr = newPtr;
				//Grow hash table if needed
				rebuildHashTable(degree, capacity);
			} else { // Type 3 => Type 3
				Neigh* __restrict newPtr = (Neigh*)globalAllocator.allocPow2(capacity * sizeof(Neigh));
				memcpy(newPtr, etype.type3.neighArr, degree * sizeof(Neigh));
				globalAllocator.freePow2(etype.type3.neighArr, degree * sizeof(Neigh));
				etype.type3.neighArr = newPtr;
				//Grow hash table if needed
				rebuildHashTable(degree, capacity);
			}
			
		}

		Neigh* __restrict currNeighArr;

		if(capacity <= TH0){
			currNeighArr = etype.type1.neigh;
		}
		else if (capacity <= TH1){
			currNeighArr = etype.type2.neighArr;
		} 
		else {
			currNeighArr = etype.type3.neighArr;
		}

		//search and insert if not found
		if(capacity <= TH1){
			//Type 1 or 2, do linear search
			for(u64 i = 0; i < degree; i++){
				if(currNeighArr[i].node == dstId){
					//found same edge, just update
					currNeighArr[i].setWeight(weight);
					return;
				}
			}
		}
		else{
			//type 3, use hash table + adj list
			u32 idx = dstId & (capacity * 2 - 1);
			DstLocPair* __restrict locMap = etype.type3.mapArr;
			DstLocPair* __restrict insLoc = nullptr;
			//probe = 0;
			while(true){
				//probe++;
				if(locMap[idx].dst == FLAG_EMPTY_SLOT){
					//edge not found, insert
					if(insLoc){
						locMap = insLoc;	//points to the first tomb stone found
					}
					locMap[idx].dst = dstId;
					locMap[idx].loc = degree;
					break;
				}
				else if((locMap[idx].dst == FLAG_TOMB_STONE) && (insLoc == nullptr)){
					insLoc = locMap + idx;
				}
				else if(locMap[idx].dst == dstId){
					//edge found, update weight
					currNeighArr[locMap[idx].loc].setWeight(weight);
					//probingDist[probe]++;
					return;
				}
				//move on
				idx++;
				if(idx == (capacity * 2)){
					idx = 0;
				}
			}
		}
		//not found, insert
		currNeighArr[degree].node = dstId;
		currNeighArr[degree].setWeight(weight);
		degree++;
		edgeCnt++;
	}


	void deleteEdge(const Idx dstId, u64& edgeCnt){
		Neigh* __restrict currNeighArr;
		Neigh* __restrict nn = nullptr;

		if(capacity <= TH0){
			currNeighArr = etype.type1.neigh;
		} 
		else if (capacity <= TH1){
			currNeighArr = etype.type2.neighArr;
		} 
		else {
			currNeighArr = etype.type3.neighArr;
		}

		//search
		if(capacity <= TH1){
			//Type 1 or 2, do linear search
			for(u64 i = 0; i < degree; i++){
				if(currNeighArr[i].node == dstId){
					nn = currNeighArr + i;
					break;
				}
			}
			if(__builtin_expect(nn != nullptr, true)){
				//edge found, delete
				degree--;
				edgeCnt--;
				nn->node = currNeighArr[degree].node;
				nn->setWeight(currNeighArr[degree].getWeight());
			}
			else{
				//edge not found, nothing to do
				return;
			}
		}
		else{
			//using hashed mode
			u32 idx = dstId & (capacity * 2 - 1);
			DstLocPair* __restrict locMap = etype.type3.mapArr;
			while(true){
				if(locMap[idx].dst == dstId){
					//edge found, delete
					degree--;
					edgeCnt--;
					//delSucc++;
					locMap[idx].dst = FLAG_TOMB_STONE; 				//invalidate previous hash-table entry

					const u32 loc = locMap[idx].loc;
					if(__builtin_expect(loc != degree, true)){		//nothing to do if last entry is removed
						const u32 node = currNeighArr[degree].node;
						//copy last entry
						currNeighArr[loc] = currNeighArr[degree];

						//point to correct location of the swapped entry
						u32 idxMoved = node & (capacity * 2 - 1);
						while(locMap[idxMoved].dst != node){
							idxMoved++;
							if(idxMoved == (capacity * 2)){
								idxMoved = 0;
							}
						}
						locMap[idxMoved].loc = loc;
					}
					break;
				}
				else if (locMap[idx].dst == FLAG_EMPTY_SLOT) {
					//edge not found, return
					return;
				}
				//move on
				idx++;
				if(idx == (capacity * 2)){
					idx = 0;
				}
			}
		}

		if((capacity > TH0) && ((degree * 4) <= capacity)){
			//time to reduce capacity
			const u64 oldCap = capacity;
			const u64 newCap = capacity / 2;
			capacity = newCap;

			Neigh* __restrict oldPtr = etype.type3.neighArr;
			Neigh* __restrict newPtr;

			if(newCap <= TH0){
				//moving from type 2 or 3 to type 1
				newPtr = etype.type1.neigh;
				capacity = TH0;
			}
			else{
				etype.type3.neighArr = (Neigh*)globalAllocator.allocPow2(newCap * sizeof(Neigh));
				newPtr = etype.type3.neighArr;
			}

			//copy old adjList and free
			memcpy(newPtr, oldPtr, degree * sizeof(Neigh));
			globalAllocator.freePow2(oldPtr, oldCap * sizeof(Neigh));

			//shrink or delete hash table if needed
			rebuildHashTable(oldCap, newCap);
		}
	}
};


template <typename Neigh>
class Vertex{
public:
	 EdgeArray<Neigh>		inEdges;
	 EdgeArray<Neigh>		outEdges;
};


#endif


#if defined(USE_GT_BALANCED_TYPE3_ONLY)

#define 	FLAG_EMPTY_SLOT			0xFFFFFFFFU
#define 	FLAG_TOMB_STONE			0xFFFFFFFEU
#define		CACHE_LINE_SIZE			64

template <typename Neigh>
class alignas(CACHE_LINE_SIZE) EdgeArray{

private:
	void rebuildHashTable(u64 oldCap, u64 newCap){
		if(mapArr){
			//free old map
			globalAllocator.freePow2(mapArr, oldCap * 2 * sizeof(DstLocPair));
			mapArr = nullptr;
		}

		//allocate new map
		mapArr = (DstLocPair*)globalAllocator.allocate(newCap * 2 * sizeof(DstLocPair));

		DstLocPair* __restrict locMap = mapArr;
		memset(locMap, -1, newCap * 2 * sizeof(DstLocPair));

		const u32 mask = newCap * 2 - 1;

		//add existing nodes to hash
		const Neigh* __restrict nn = neighArr;
		for(u64 i = 0; i < degree; i++){
			const u32 dst = nn[i].node;
			u32 idx = dst & mask;
			while(true){
				if(locMap[idx].dst == FLAG_EMPTY_SLOT){
					//found insertion point
					locMap[idx].dst = dst;
					locMap[idx].loc = i;
					break;
				}
				//move on
				idx++;
				if(idx == (newCap * 2)){
					idx = 0;
				}
			}
		}

	}


public:

	u64 degree = 0;
	u64 capacity = 0;

	Neigh* 			__restrict neighArr = nullptr;
	DstLocPair* 	__restrict mapArr = nullptr;

	void insertEdge(const Idx dstId, const Weight weight, u64& edgeCnt){
		//First, check if needs expanding
		if(__builtin_expect(degree == capacity, false)){
			capacity = getNextPow2MinRet(capacity * 2);
			Neigh* __restrict newPtr = (Neigh*)globalAllocator.allocPow2(capacity * sizeof(Neigh));

			if(neighArr){
				memcpy(newPtr, neighArr, degree * sizeof(Neigh));
				globalAllocator.freePow2(neighArr, capacity / 2 * sizeof(Neigh));
			}
			neighArr = newPtr;

			//Grow hash table if needed
			rebuildHashTable(capacity / 2, capacity);
		}

		//type 3, use hash table + adj list
		u32 idx = dstId & (capacity * 2 - 1);
		DstLocPair* __restrict locMap = mapArr;
		DstLocPair* __restrict insLoc = nullptr;
		//probe = 0;
		while(true){
			//probe++;
			if(locMap[idx].dst == FLAG_EMPTY_SLOT){
				//edge not found, insert
				if(insLoc){
					locMap = insLoc;	//points to the first tomb stone found
				}
				locMap[idx].dst = dstId;
				locMap[idx].loc = degree;
				break;
			}
			else if((locMap[idx].dst == FLAG_TOMB_STONE) && (insLoc == nullptr)){
				insLoc = locMap + idx;
			}
			else if(locMap[idx].dst == dstId){
				//edge found, update weight
				neighArr[locMap[idx].loc].setWeight(weight);
				//probingDist[probe]++;
				return;
			}
			//move on
			idx++;
			if(idx == (capacity * 2)){
				idx = 0;
			}
		}

		//not found, insert
		neighArr[degree].node = dstId;
		neighArr[degree].setWeight(weight);
		degree++;
		edgeCnt++;
	}


	void deleteEdge(const Idx dstId, u64& edgeCnt){
		Neigh* __restrict currNeighArr;
		Neigh* __restrict nn = nullptr;

		//using hashed mode
		u32 idx = dstId & (capacity * 2 - 1);
		DstLocPair* __restrict locMap = mapArr;
		while(true){
			if(locMap[idx].dst == dstId){
				//edge found, delete
				degree--;
				edgeCnt--;
				//delSucc++;
				locMap[idx].dst = FLAG_TOMB_STONE; 				//invalidate previous hash-table entry

				const u32 loc = locMap[idx].loc;
				if(__builtin_expect(loc != degree, true)){		//nothing to do if last entry is removed
					const u32 node = currNeighArr[degree].node;
					//copy last entry
					currNeighArr[loc] = currNeighArr[degree];

					//point to correct location of the swapped entry
					u32 idxMoved = node & (capacity * 2 - 1);
					while(locMap[idxMoved].dst != node){
						idxMoved++;
						if(idxMoved == (capacity * 2)){
							idxMoved = 0;
						}
					}
					locMap[idxMoved].loc = loc;
				}
				break;
			}
			else if (locMap[idx].dst == FLAG_EMPTY_SLOT) {
				//edge not found, return
				return;
			}
			//move on
			idx++;
			if(idx == (capacity * 2)){
				idx = 0;
			}
		}

		if(degree * 4 <= capacity){
			//time to reduce capacity
			const u64 oldCap = capacity;
			const u64 newCap = capacity / 2;
			capacity = newCap;

			Neigh* __restrict oldPtr = neighArr;
			Neigh* __restrict newPtr = (Neigh*)globalAllocator.allocPow2(newCap * sizeof(Neigh));

			neighArr = newPtr;

			//copy old adjList and free
			memcpy(newPtr, oldPtr, degree * sizeof(Neigh));
			globalAllocator.freePow2(oldPtr, oldCap * sizeof(Neigh));

			//shrink or delete hash table if needed
			rebuildHashTable(oldCap, newCap);
		}
	}
};


template <typename Neigh>
class Vertex{
public:
	 EdgeArray<Neigh>		inEdges;
	 EdgeArray<Neigh>		outEdges;
};


#endif


#ifdef USE_GT_BALANCED_STDMAP

#define		CACHE_LINE_SIZE			64

template <typename Neigh>
class alignas(CACHE_LINE_SIZE) EdgeArray{
public:

	const static u64 TH0 = ((CACHE_LINE_SIZE - sizeof(u32) - sizeof(u32)) / sizeof(Neigh));
	const static u64 TH1 = HYBRID_HASH_PARTITION;

	u32 degree = 0;
	u32 capacity = TH0;

	union {
		struct {
			Neigh neigh[TH0];
		} type1;

		struct {
			Neigh* 			__restrict neighArr = nullptr;
			graphite_hashmap* 	__restrict locMap = nullptr;
		} type2_3;
	} etype;

	u8 __pad[CACHE_LINE_SIZE - sizeof(u32) - sizeof(u32) - sizeof(etype)];

	void insertEdge(const Idx dstId, const Weight weight, u64& edgeCnt){
		//First, check if needs expanding
		if(__builtin_expect(degree == capacity, false)){
			capacity = getNextPow2(capacity * 2);
			Neigh* __restrict newPtr = (Neigh*)globalAllocator.allocPow2(capacity * sizeof(Neigh));

			if(degree <= TH0){	//Going from Type 1 to Type 2
				memcpy(newPtr, etype.type1.neigh, degree * sizeof(Neigh));
				etype.type2_3.locMap = nullptr;
			}
			else{				//Type 2 or 3
				memcpy(newPtr, etype.type2_3.neighArr, degree * sizeof(Neigh));
				globalAllocator.freePow2(etype.type2_3.neighArr, capacity / 2 * sizeof(Neigh));
			}
			etype.type2_3.neighArr = newPtr;

			if((capacity > TH1) && (etype.type2_3.locMap == nullptr)){
				//type 3 and locMap not yet allocated
				etype.type2_3.locMap = (graphite_hashmap*)globalAllocator.allocate(sizeof(graphite_hashmap));
				new (etype.type2_3.locMap) graphite_hashmap();	//placement new
			}
		}

		Neigh* __restrict currNeighArr;

		if(capacity <= TH0){
			currNeighArr = etype.type1.neigh;
		}
		else{
			currNeighArr = etype.type2_3.neighArr;
		}

		//search and insert if not found
		if(capacity <= TH1){
			//Type 1 or 2, do linear search
			for(u64 i = 0; i < degree; i++){
				if(currNeighArr[i].node == dstId){
					//found same edge, just update
					currNeighArr[i].setWeight(weight);
					return;
				}
			}
		}
		else{
			//type 3, use hash table + adj list
			graphite_hashmap* __restrict locMap = etype.type2_3.locMap;

			const auto& iter = locMap->find(dstId);
			if(iter != locMap->end()){
				//found same edge, just update
				currNeighArr[iter->second].setWeight(weight);
				return;
			}

			//edge not found, must insert
			locMap->insert({dstId, degree});
		}
		//not found, insert
		currNeighArr[degree].node = dstId;
		currNeighArr[degree].setWeight(weight);
		degree++;
		edgeCnt++;
	}


	void deleteEdge(const Idx dstId, u64& edgeCnt){
		Neigh* __restrict currNeighArr;
		Neigh* __restrict nn = nullptr;

		if(capacity <= TH0){
			currNeighArr = etype.type1.neigh;
		}
		else{
			currNeighArr = etype.type2_3.neighArr;
		}

		//search
		if(capacity <= TH1){
			//Type 1 or 2, do linear search
			for(u64 i = 0; i < degree; i++){
				if(currNeighArr[i].node == dstId){
					nn = currNeighArr + i;
					break;
				}
			}
			if(__builtin_expect(nn != nullptr, true)){
				//edge found, delete
				degree--;
				edgeCnt--;
				nn->node = currNeighArr[degree].node;
				nn->setWeight(currNeighArr[degree].getWeight());
			}
			else{
				//edge not found, nothing to do
				return;
			}
		}
		else{
			//using hashed mode
			graphite_hashmap* __restrict locMap = etype.type2_3.locMap;
			const auto& iter = locMap->find(dstId);
			if(iter == locMap->end()){
				//not found, nothing to do
				return;
			}

			//edge found
			degree--;
			edgeCnt--;
			const u32 loc = iter->second;
			locMap->erase(iter);
			if(__builtin_expect(loc != degree, true)){
				const u32 node = currNeighArr[degree].node;
				//copy last entry
				currNeighArr[loc] = currNeighArr[degree];
				const auto &nodeIter = locMap->find(node);
				nodeIter->second = loc;
			}
		}

		if((capacity > TH0) && ((degree * 4) <= capacity)){
			//time to reduce capacity
			const u64 oldCap = capacity;
			const u64 newCap = capacity / 2;
			capacity = newCap;

			Neigh* __restrict oldPtr = etype.type2_3.neighArr;
			Neigh* __restrict newPtr;

			if(newCap <= TH0){
				//moving from type 2 or 3 to type 1
				newPtr = etype.type1.neigh;
				capacity = TH0;
			}
			else{
				etype.type2_3.neighArr = (Neigh*)globalAllocator.allocPow2(newCap * sizeof(Neigh));
				newPtr = etype.type2_3.neighArr;
			}

			//copy old adjList and free
			memcpy(newPtr, oldPtr, degree * sizeof(Neigh));
			globalAllocator.freePow2(oldPtr, oldCap * sizeof(Neigh));

			if((oldCap > TH1) && (newCap <= TH1)){
				//moving from type 3 to type 2/1
				globalAllocator.freePow2(etype.type2_3.locMap, getNextPow2(sizeof(graphite_hashmap)));
				etype.type2_3.locMap = nullptr;
			}
		}
	}
};


template <typename Neigh>
class Vertex{
public:
	 EdgeArray<Neigh>		inEdges;
	 EdgeArray<Neigh>		outEdges;
};


#endif



#if defined(USE_GT_BALANCED_MALLOC_STDMAP) || defined(USE_GT_BALANCED_RHH)

#define		CACHE_LINE_SIZE			64

template <typename Neigh>
class alignas(CACHE_LINE_SIZE) EdgeArray{
public:

	const static u64 TH0 = ((CACHE_LINE_SIZE - sizeof(u32) - sizeof(u32)) / sizeof(Neigh));
	const static u64 TH1 = HYBRID_HASH_PARTITION;

	u32 degree = 0;
	u32 capacity = TH0;

	union {
		struct {
			Neigh neigh[TH0];
		} type1;

		struct {
			Neigh* 			__restrict neighArr = nullptr;
			graphite_hashmap* 	__restrict locMap = nullptr;
		} type2_3;
	} etype;

	u8 __pad[CACHE_LINE_SIZE - sizeof(u32) - sizeof(u32) - sizeof(etype)];

	void insertEdge(const Idx dstId, const Weight weight, u64& edgeCnt){
		//First, check if needs expanding
		if(__builtin_expect(degree == capacity, false)){
			capacity = getNextPow2(capacity * 2);
			Neigh* __restrict newPtr = (Neigh*)aligned_alloc(64, capacity * sizeof(Neigh));

			if(degree <= TH0){	//Going from Type 1 to Type 2
				memcpy(newPtr, etype.type1.neigh, degree * sizeof(Neigh));
				etype.type2_3.locMap = nullptr;
			}
			else{				//Type 2 or 3
				memcpy(newPtr, etype.type2_3.neighArr, degree * sizeof(Neigh));
				free(etype.type2_3.neighArr);
			}
			etype.type2_3.neighArr = newPtr;

			if((capacity > TH1) && (etype.type2_3.locMap == nullptr)){
				//type 3 and locMap not yet allocated
				etype.type2_3.locMap = (graphite_hashmap*)malloc(sizeof(graphite_hashmap));
				new (etype.type2_3.locMap) graphite_hashmap();	//placement new
			}
		}

		Neigh* __restrict currNeighArr;

		if(capacity <= TH0){
			currNeighArr = etype.type1.neigh;
		}
		else{
			currNeighArr = etype.type2_3.neighArr;
		}

		//search and insert if not found
		if(capacity <= TH1){
			//Type 1 or 2, do linear search
			for(u64 i = 0; i < degree; i++){
				if(currNeighArr[i].node == dstId){
					//found same edge, just update
					currNeighArr[i].setWeight(weight);
					return;
				}
			}
		}
		else{
			//type 3, use hash table + adj list
			graphite_hashmap* __restrict locMap = etype.type2_3.locMap;

			const auto& iter = locMap->find(dstId);
			if(iter != locMap->end()){
				//found same edge, just update
				currNeighArr[iter->second].setWeight(weight);
				return;
			}

			//edge not found, must insert
			locMap->insert({dstId, degree});
		}
		//not found, insert
		currNeighArr[degree].node = dstId;
		currNeighArr[degree].setWeight(weight);
		degree++;
		edgeCnt++;
	}


	void deleteEdge(const Idx dstId, u64& edgeCnt){
		Neigh* __restrict currNeighArr;
		Neigh* __restrict nn = nullptr;

		if(capacity <= TH0){
			currNeighArr = etype.type1.neigh;
		}
		else{
			currNeighArr = etype.type2_3.neighArr;
		}

		//search
		if(capacity <= TH1){
			//Type 1 or 2, do linear search
			for(u64 i = 0; i < degree; i++){
				if(currNeighArr[i].node == dstId){
					nn = currNeighArr + i;
					break;
				}
			}
			if(__builtin_expect(nn != nullptr, true)){
				//edge found, delete
				degree--;
				edgeCnt--;
				nn->node = currNeighArr[degree].node;
				nn->setWeight(currNeighArr[degree].getWeight());
			}
			else{
				//edge not found, nothing to do
				return;
			}
		}
		else{
			//using hashed mode
			graphite_hashmap* __restrict locMap = etype.type2_3.locMap;
			const auto& iter = locMap->find(dstId);
			if(iter == locMap->end()){
				//not found, nothing to do
				return;
			}

			//edge found
			degree--;
			edgeCnt--;
			const u32 loc = iter->second;
			locMap->erase(iter);
			if(__builtin_expect(loc != degree, true)){
				const u32 node = currNeighArr[degree].node;
				//copy last entry
				currNeighArr[loc] = currNeighArr[degree];
				const auto &nodeIter = locMap->find(node);
				nodeIter->second = loc;
			}
		}

		if((capacity > TH0) && ((degree * 4) <= capacity)){
			//time to reduce capacity
			const u64 oldCap = capacity;
			const u64 newCap = capacity / 2;
			capacity = newCap;

			Neigh* __restrict oldPtr = etype.type2_3.neighArr;
			Neigh* __restrict newPtr;

			if(newCap <= TH0){
				//moving from type 2 or 3 to type 1
				newPtr = etype.type1.neigh;
				capacity = TH0;
			}
			else{
				etype.type2_3.neighArr = (Neigh*)aligned_alloc(64, newCap * sizeof(Neigh));
				newPtr = etype.type2_3.neighArr;
			}

			//copy old adjList and free
			memcpy(newPtr, oldPtr, degree * sizeof(Neigh));
			free(oldPtr);

			if((oldCap > TH1) && (newCap <= TH1)){
				//moving from type 3 to type 2/1
				free(etype.type2_3.locMap);
				etype.type2_3.locMap = nullptr;
			}
		}
	}
};


template <typename Neigh>
class Vertex{
public:
	 EdgeArray<Neigh>		inEdges;
	 EdgeArray<Neigh>		outEdges;
};


#endif



#if defined(USE_GT_BALANCED_TSL_RHH)

#define		CACHE_LINE_SIZE			64

template <typename Neigh>
class alignas(CACHE_LINE_SIZE) EdgeArray{
public:

	const static u64 TH0 = ((CACHE_LINE_SIZE - sizeof(u32) - sizeof(u32)) / sizeof(Neigh));
	const static u64 TH1 = HYBRID_HASH_PARTITION;

	u32 degree = 0;
	u32 capacity = TH0;

	union {
		struct {
			Neigh neigh[TH0];
		} type1;

		struct {
			Neigh* 			__restrict neighArr = nullptr;
			graphite_hashmap* 	__restrict locMap = nullptr;
		} type2_3;
	} etype;

	u8 __pad[CACHE_LINE_SIZE - sizeof(u32) - sizeof(u32) - sizeof(etype)];

	void insertEdge(const Idx dstId, const Weight weight, u64& edgeCnt){
		//First, check if needs expanding
		if(__builtin_expect(degree == capacity, false)){
			capacity = getNextPow2(capacity * 2);
			Neigh* __restrict newPtr = (Neigh*)aligned_alloc(64, capacity * sizeof(Neigh));

			if(degree <= TH0){	//Going from Type 1 to Type 2
				memcpy(newPtr, etype.type1.neigh, degree * sizeof(Neigh));
				etype.type2_3.locMap = nullptr;
			}
			else{				//Type 2 or 3
				memcpy(newPtr, etype.type2_3.neighArr, degree * sizeof(Neigh));
				free(etype.type2_3.neighArr);
			}
			etype.type2_3.neighArr = newPtr;

			if((capacity > TH1) && (etype.type2_3.locMap == nullptr)){
				//type 3 and locMap not yet allocated
				etype.type2_3.locMap = (graphite_hashmap*)malloc(sizeof(graphite_hashmap));
				new (etype.type2_3.locMap) graphite_hashmap();	//placement new
			}
		}

		Neigh* __restrict currNeighArr;

		if(capacity <= TH0){
			currNeighArr = etype.type1.neigh;
		}
		else{
			currNeighArr = etype.type2_3.neighArr;
		}

		//search and insert if not found
		if(capacity <= TH1){
			//Type 1 or 2, do linear search
			for(u64 i = 0; i < degree; i++){
				if(currNeighArr[i].node == dstId){
					//found same edge, just update
					currNeighArr[i].setWeight(weight);
					return;
				}
			}
		}
		else{
			//type 3, use hash table + adj list
			graphite_hashmap* __restrict locMap = etype.type2_3.locMap;

			const auto& iter = locMap->find(dstId);
			if(iter != locMap->end()){
				//found same edge, just update
				currNeighArr[iter->second].setWeight(weight);
				return;
			}

			//edge not found, must insert
			locMap->insert({dstId, degree});
		}
		//not found, insert
		currNeighArr[degree].node = dstId;
		currNeighArr[degree].setWeight(weight);
		degree++;
		edgeCnt++;
	}


	void deleteEdge(const Idx dstId, u64& edgeCnt){
		Neigh* __restrict currNeighArr;
		Neigh* __restrict nn = nullptr;

		if(capacity <= TH0){
			currNeighArr = etype.type1.neigh;
		}
		else{
			currNeighArr = etype.type2_3.neighArr;
		}

		//search
		if(capacity <= TH1){
			//Type 1 or 2, do linear search
			for(u64 i = 0; i < degree; i++){
				if(currNeighArr[i].node == dstId){
					nn = currNeighArr + i;
					break;
				}
			}
			if(__builtin_expect(nn != nullptr, true)){
				//edge found, delete
				degree--;
				edgeCnt--;
				nn->node = currNeighArr[degree].node;
				nn->setWeight(currNeighArr[degree].getWeight());
			}
			else{
				//edge not found, nothing to do
				return;
			}
		}
		else{
			//using hashed mode
			graphite_hashmap* __restrict locMap = etype.type2_3.locMap;
			const auto& iter = locMap->find(dstId);
			if(iter == locMap->end()){
				//not found, nothing to do
				return;
			}

			//edge found
			degree--;
			edgeCnt--;
			const u32 loc = iter->second;
			locMap->erase(iter);
			if(__builtin_expect(loc != degree, true)){
				const u32 node = currNeighArr[degree].node;
				//copy last entry
				currNeighArr[loc] = currNeighArr[degree];
				(*locMap)[node] = loc;
				//const auto &nodeIter = locMap->find(node);
				//nodeIter->second = loc;
			}
		}

		if((capacity > TH0) && ((degree * 4) <= capacity)){
			//time to reduce capacity
			const u64 oldCap = capacity;
			const u64 newCap = capacity / 2;
			capacity = newCap;

			Neigh* __restrict oldPtr = etype.type2_3.neighArr;
			Neigh* __restrict newPtr;

			if(newCap <= TH0){
				//moving from type 2 or 3 to type 1
				newPtr = etype.type1.neigh;
				capacity = TH0;
			}
			else{
				etype.type2_3.neighArr = (Neigh*)aligned_alloc(64, newCap * sizeof(Neigh));
				newPtr = etype.type2_3.neighArr;
			}

			//copy old adjList and free
			memcpy(newPtr, oldPtr, degree * sizeof(Neigh));
			free(oldPtr);

			if((oldCap > TH1) && (newCap <= TH1)){
				//moving from type 3 to type 2/1
				free(etype.type2_3.locMap);
				etype.type2_3.locMap = nullptr;
			}
		}
	}
};


template <typename Neigh>
class Vertex{
public:
	 EdgeArray<Neigh>		inEdges;
	 EdgeArray<Neigh>		outEdges;
};


#endif




#ifdef USE_GT_BALANCED_ABSEIL

#define		CACHE_LINE_SIZE			64

template <typename Neigh>
class alignas(CACHE_LINE_SIZE) EdgeArray{
public:

	const static u64 TH0 = ((CACHE_LINE_SIZE - sizeof(u32) - sizeof(u32)) / sizeof(Neigh));
	const static u64 TH1 = HYBRID_HASH_PARTITION;

	u32 degree = 0;
	u32 capacity = TH0;

	union {
		struct {
			Neigh neigh[TH0];
		} type1;

		struct {
			Neigh* 			__restrict neighArr = nullptr;
			graphite_hashmap* 	__restrict locMap = nullptr;
		} type2_3;
	} etype;

	u8 __pad[CACHE_LINE_SIZE - sizeof(u32) - sizeof(u32) - sizeof(etype)];

	void insertEdge(const Idx dstId, const Weight weight, u64& edgeCnt){
		//First, check if needs expanding
		if(__builtin_expect(degree == capacity, false)){
			capacity = getNextPow2(capacity * 2);
			Neigh* __restrict newPtr = (Neigh*)aligned_alloc(64, capacity * sizeof(Neigh));

			if(degree <= TH0){	//Going from Type 1 to Type 2
				memcpy(newPtr, etype.type1.neigh, degree * sizeof(Neigh));
				etype.type2_3.locMap = nullptr;
			}
			else{				//Type 2 or 3
				memcpy(newPtr, etype.type2_3.neighArr, degree * sizeof(Neigh));
				free(etype.type2_3.neighArr);
			}
			etype.type2_3.neighArr = newPtr;

			if((capacity > TH1) && (etype.type2_3.locMap == nullptr)){
				//type 3 and locMap not yet allocated
				etype.type2_3.locMap = (graphite_hashmap*)malloc(sizeof(graphite_hashmap));
				new (etype.type2_3.locMap) graphite_hashmap();	//placement new
			}
		}

		Neigh* __restrict currNeighArr;

		if(capacity <= TH0){
			currNeighArr = etype.type1.neigh;
		}
		else{
			currNeighArr = etype.type2_3.neighArr;
		}

		//search and insert if not found
		if(capacity <= TH1){
			//Type 1 or 2, do linear search
			for(u64 i = 0; i < degree; i++){
				if(currNeighArr[i].node == dstId){
					//found same edge, just update
					currNeighArr[i].setWeight(weight);
					return;
				}
			}
		}
		else{
			//type 3, use hash table + adj list
			graphite_hashmap* __restrict locMap = etype.type2_3.locMap;

			const auto& iter = locMap->find(dstId);
			if(iter != locMap->end()){
				//found same edge, just update
				currNeighArr[iter->second].setWeight(weight);
				return;
			}

			//edge not found, must insert
			locMap->insert({dstId, degree});
		}
		//not found, insert
		currNeighArr[degree].node = dstId;
		currNeighArr[degree].setWeight(weight);
		degree++;
		edgeCnt++;
	}


	void deleteEdge(const Idx dstId, u64& edgeCnt){
		Neigh* __restrict currNeighArr;
		Neigh* __restrict nn = nullptr;

		if(capacity <= TH0){
			currNeighArr = etype.type1.neigh;
		}
		else{
			currNeighArr = etype.type2_3.neighArr;
		}

		//search
		if(capacity <= TH1){
			//Type 1 or 2, do linear search
			for(u64 i = 0; i < degree; i++){
				if(currNeighArr[i].node == dstId){
					nn = currNeighArr + i;
					break;
				}
			}
			if(__builtin_expect(nn != nullptr, true)){
				//edge found, delete
				degree--;
				edgeCnt--;
				nn->node = currNeighArr[degree].node;
				nn->setWeight(currNeighArr[degree].getWeight());
			}
			else{
				//edge not found, nothing to do
				return;
			}
		}
		else{
			//using hashed mode
			graphite_hashmap* __restrict locMap = etype.type2_3.locMap;
			const auto& iter = locMap->find(dstId);
			if(iter == locMap->end()){
				//not found, nothing to do
				return;
			}

			//edge found
			degree--;
			edgeCnt--;
			const u32 loc = iter->second;
			locMap->erase(iter);
			if(__builtin_expect(loc != degree, true)){
				const u32 node = currNeighArr[degree].node;
				//copy last entry
				currNeighArr[loc] = currNeighArr[degree];
				const auto &nodeIter = locMap->find(node);
				nodeIter->second = loc;
			}
		}

		if((capacity > TH0) && ((degree * 4) <= capacity)){
			//time to reduce capacity
			const u64 oldCap = capacity;
			const u64 newCap = capacity / 2;
			capacity = newCap;

			Neigh* __restrict oldPtr = etype.type2_3.neighArr;
			Neigh* __restrict newPtr;

			if(newCap <= TH0){
				//moving from type 2 or 3 to type 1
				newPtr = etype.type1.neigh;
				capacity = TH0;
			}
			else{
				etype.type2_3.neighArr = (Neigh*)aligned_alloc(64, newCap * sizeof(Neigh));
				newPtr = etype.type2_3.neighArr;
			}

			//copy old adjList and free
			memcpy(newPtr, oldPtr, degree * sizeof(Neigh));
			free(oldPtr);

			if((oldCap > TH1) && (newCap <= TH1)){
				//moving from type 3 to type 2/1
				free(etype.type2_3.locMap);
				etype.type2_3.locMap = nullptr;
			}
		}
	}
};


template <typename Neigh>
class Vertex{
public:
	 EdgeArray<Neigh>		inEdges;
	 EdgeArray<Neigh>		outEdges;
};


#endif



#ifdef USE_GT_BALANCED_MALLOC

#define 	FLAG_EMPTY_SLOT			0xFFFFFFFFU
#define 	FLAG_TOMB_STONE			0xFFFFFFFEU
#define		CACHE_LINE_SIZE			64

template <typename Neigh>
class alignas(CACHE_LINE_SIZE) EdgeArray{

private:
	void rebuildHashTable(u64 oldCap, u64 newCap){
		if(oldCap > TH1){
			//free old map
			free(etype.type2_3.mapArr);
			etype.type2_3.mapArr = nullptr;
		}

		if(newCap > TH1){
			//allocate new map
			etype.type2_3.mapArr = (DstLocPair*)aligned_alloc(64, newCap * 2 * sizeof(DstLocPair));

			DstLocPair* __restrict locMap = etype.type2_3.mapArr;
			memset(locMap, -1, newCap * 2 * sizeof(DstLocPair));

			const u32 mask = newCap * 2 - 1;

			//add existing nodes to hash
			const Neigh* __restrict nn = etype.type2_3.neighArr;
			for(u64 i = 0; i < degree; i++){
				const u32 dst = nn[i].node;
				u32 idx = dst & mask;
				while(true){
					if(locMap[idx].dst == FLAG_EMPTY_SLOT){
						//found insertion point
						locMap[idx].dst = dst;
						locMap[idx].loc = i;
						break;
					}
					//move on
					idx++;
					if(idx == (newCap * 2)){
						idx = 0;
					}
				}
			}
		}
	}


public:

	const static u64 TH0 = ((CACHE_LINE_SIZE - sizeof(u32) - sizeof(u32)) / sizeof(Neigh));
	const static u64 TH1 = HYBRID_HASH_PARTITION;

	u32 degree = 0;
	u32 capacity = TH0;

	union {
		struct {
			Neigh neigh[TH0];
		} type1;

		struct {
			Neigh* 			__restrict neighArr = nullptr;
			DstLocPair* 	__restrict mapArr = nullptr;
		} type2_3;
	} etype;

	u8 __pad[CACHE_LINE_SIZE - sizeof(u32) - sizeof(u32) - sizeof(etype)];

	void insertEdge(const Idx dstId, const Weight weight, u64& edgeCnt){
		//First, check if needs expanding
		if(__builtin_expect(degree == capacity, false)){
			capacity = getNextPow2(capacity * 2);
			Neigh* __restrict newPtr = (Neigh*)aligned_alloc(64, capacity * sizeof(Neigh));

			if(degree <= TH0){	//Going from Type 1 to Type 2
				memcpy(newPtr, etype.type1.neigh, degree * sizeof(Neigh));
				etype.type2_3.mapArr = nullptr;
			}
			else{				//Type 2 or 3
				memcpy(newPtr, etype.type2_3.neighArr, degree * sizeof(Neigh));
				free(etype.type2_3.neighArr);
			}
			etype.type2_3.neighArr = newPtr;

			//Grow hash table if needed
			rebuildHashTable(capacity / 2, capacity);
		}

		Neigh* __restrict currNeighArr;

		if(capacity <= TH0){
			currNeighArr = etype.type1.neigh;
		}
		else{
			currNeighArr = etype.type2_3.neighArr;
		}

		//search and insert if not found
		if(capacity <= TH1){
			//Type 1 or 2, do linear search
			for(u64 i = 0; i < degree; i++){
				if(currNeighArr[i].node == dstId){
					//found same edge, just update
					currNeighArr[i].setWeight(weight);
					return;
				}
			}
		}
		else{
			//type 3, use hash table + adj list
			u32 idx = dstId & (capacity * 2 - 1);
			DstLocPair* __restrict locMap = etype.type2_3.mapArr;
			DstLocPair* __restrict insLoc = nullptr;
			//probe = 0;
			while(true){
				//probe++;
				if(locMap[idx].dst == FLAG_EMPTY_SLOT){
					//edge not found, insert
					if(insLoc){
						locMap = insLoc;	//points to the first tomb stone found
					}
					locMap[idx].dst = dstId;
					locMap[idx].loc = degree;
					break;
				}
				else if((locMap[idx].dst == FLAG_TOMB_STONE) && (insLoc == nullptr)){
					insLoc = locMap + idx;
				}
				else if(locMap[idx].dst == dstId){
					//edge found, update weight
					currNeighArr[locMap[idx].loc].setWeight(weight);
					//probingDist[probe]++;
					return;
				}
				//move on
				idx++;
				if(idx == (capacity * 2)){
					idx = 0;
				}
			}
		}
		//not found, insert
		currNeighArr[degree].node = dstId;
		currNeighArr[degree].setWeight(weight);
		degree++;
		edgeCnt++;
	}


	void deleteEdge(const Idx dstId, u64& edgeCnt){
		Neigh* __restrict currNeighArr;
		Neigh* __restrict nn = nullptr;

		if(capacity <= TH0){
			currNeighArr = etype.type1.neigh;
		}
		else{
			currNeighArr = etype.type2_3.neighArr;
		}

		//search
		if(capacity <= TH1){
			//Type 1 or 2, do linear search
			for(u64 i = 0; i < degree; i++){
				if(currNeighArr[i].node == dstId){
					nn = currNeighArr + i;
					break;
				}
			}
			if(__builtin_expect(nn != nullptr, true)){
				//edge found, delete
				degree--;
				edgeCnt--;
				nn->node = currNeighArr[degree].node;
				nn->setWeight(currNeighArr[degree].getWeight());
			}
			else{
				//edge not found, nothing to do
				return;
			}
		}
		else{
			//using hashed mode
			u32 idx = dstId & (capacity * 2 - 1);
			DstLocPair* __restrict locMap = etype.type2_3.mapArr;
			while(true){
				if(locMap[idx].dst == dstId){
					//edge found, delete
					degree--;
					edgeCnt--;
					//delSucc++;
					locMap[idx].dst = FLAG_TOMB_STONE; 				//invalidate previous hash-table entry

					const u32 loc = locMap[idx].loc;
					if(__builtin_expect(loc != degree, true)){		//nothing to do if last entry is removed
						const u32 node = currNeighArr[degree].node;
						//copy last entry
						currNeighArr[loc] = currNeighArr[degree];

						//point to correct location of the swapped entry
						u32 idxMoved = node & (capacity * 2 - 1);
						while(locMap[idxMoved].dst != node){
							idxMoved++;
							if(idxMoved == (capacity * 2)){
								idxMoved = 0;
							}
						}
						locMap[idxMoved].loc = loc;
					}
					break;
				}
				else if (locMap[idx].dst == FLAG_EMPTY_SLOT) {
					//edge not found, return
					return;
				}
				//move on
				idx++;
				if(idx == (capacity * 2)){
					idx = 0;
				}
			}
		}

		if((capacity > TH0) && ((degree * 4) <= capacity)){
			//time to reduce capacity
			const u64 oldCap = capacity;
			const u64 newCap = capacity / 2;
			capacity = newCap;

			Neigh* __restrict oldPtr = etype.type2_3.neighArr;
			Neigh* __restrict newPtr;

			if(newCap <= TH0){
				//moving from type 2 or 3 to type 1
				newPtr = etype.type1.neigh;
				capacity = TH0;
			}
			else{
				etype.type2_3.neighArr = (Neigh*)aligned_alloc(64, newCap * sizeof(Neigh));
				newPtr = etype.type2_3.neighArr;
			}

			//copy old adjList and free
			memcpy(newPtr, oldPtr, degree * sizeof(Neigh));
			free(oldPtr);

			//shrink or delete hash table if needed
			rebuildHashTable(oldCap, newCap);
		}
	}
};


template <typename Neigh>
class Vertex{
public:
	 EdgeArray<Neigh>		inEdges;
	 EdgeArray<Neigh>		outEdges;
};


#endif



#ifdef USE_GT_UPDATE

#define 	FLAG_EMPTY_SLOT			0xFFFFFFFFU
#define 	FLAG_TOMB_STONE			0xFFFFFFFEU
#define		CACHE_LINE_SIZE			64

template <typename Neigh>
class alignas(CACHE_LINE_SIZE) EdgeArray{

private:
	void rebuildHashTable(u64 oldCap, u64 newCap){
		if(oldCap > TH1){
			//free old map
			globalAllocator.freePow2(etype.type2_3.mapArr, oldCap * 2 * sizeof(Neigh));
			etype.type2_3.mapArr = nullptr;
		}

		if(newCap > TH1){
			//allocate new map
			etype.type2_3.mapArr = (Neigh*)globalAllocator.allocate(newCap * 2 * sizeof(Neigh));

			Neigh* __restrict locMap = etype.type2_3.mapArr;
			memset(locMap, -1, newCap * 2 * sizeof(DstLocPair));

			const u32 mask = newCap * 2 - 1;

			//add existing nodes to hash
			const Neigh* __restrict nn = etype.type2_3.neighArr;
			for(u64 i = 0; i < degree; i++){
				const u32 dst = nn[i].node;
				u32 idx = dst & mask;
				while(true){
					if(locMap[idx].dst == FLAG_EMPTY_SLOT){
						//found insertion point
						locMap[idx].dst = dst;
						locMap[idx].loc = i;
						break;
					}
					//move on
					idx++;
					if(idx == (newCap * 2)){
						idx = 0;
					}
				}
			}
		}
	}


public:

	const static u64 TH0 = ((CACHE_LINE_SIZE - sizeof(u32) - sizeof(u32)) / sizeof(Neigh));
	const static u64 TH1 = HYBRID_HASH_PARTITION;

	u32 degree = 0;
	u32 capacity = TH0;

	union {
		struct {
			Neigh neigh[TH0];
		} type1;

		struct {
			Neigh* 			__restrict neighArr = nullptr;
			Neigh* 			__restrict mapArr = nullptr;
		} type2_3;
	} etype;

	u8 __pad[CACHE_LINE_SIZE - sizeof(u32) - sizeof(u32) - sizeof(etype)];

	void insertEdge(const Idx dstId, const Weight weight, u64& edgeCnt){
		//First, check if needs expanding
		if(__builtin_expect(degree == capacity, false)){
			Neigh* __restrict oldArrPtr;
			Neigh* __restrict newArrPtr;

			Neigh* __restrict oldMapPtr;
			Neigh* __restrict newMapPtr;

			const u64 oldCap = capacity;
			const u64 newCap = getNextPow2(capacity * 2);
			capacity = newCap;

			if(oldCap <= TH0){
				oldArrPtr = etype.type1.neigh;
			}
			else{
				oldArrPtr = etype.type2_3.neighArr;
			}

			if(newCap <= TH1){
				//Type 1->2 or 2->2
				//create newArrPtr and copy
				newArrPtr = (Neigh*)globalAllocator.allocPow2(newCap * sizeof(Neigh));
				memcpy(newArrPtr, oldArrPtr, degree * sizeof(Neigh));
			}
			else{
				//1->3 or 2->3 or 3->3
				//create newMapPtr
				newMapPtr = (Neigh*)globalAllocator.allocPow2(newCap * 2 * sizeof(Neigh));
				memset(newMapPtr, -1, newCap * 2 * sizeof(Neigh));

				const u32 mask = newCap * 2 - 1;

				//add existing nodes to the hash table
				if(oldArrPtr){
					//1->3 or 2->3 (or 3->3 with valid adj list)
					for(u64 i = 0; i < degree; i++){
						const u32 dst = oldArrPtr[i].node;
						u32 idx = dst & mask;
						while(true){
							if(newMapPtr[idx].node == FLAG_EMPTY_SLOT){
								//found insertion point
								newMapPtr[idx] = oldArrPtr[i];
								break;
							}
							//move on
							idx++;
							if(idx == (newCap * 2)){
								idx = 0;
							}
						}
					}

				}
				else{
					oldMapPtr = etype.type2_3.mapArr;
					//3->3 without valid adj list
					for(u64 i = 0; i < oldCap * 2; i++){
						if((u32)(oldMapPtr[i].node) < FLAG_TOMB_STONE){
							//valid edge, insert to new map
							const u32 dst = oldMapPtr[i].node;
							u32 idx = dst & mask;
							while(true){
								if(newMapPtr[idx].node == FLAG_EMPTY_SLOT){
									//found insertion point
									newMapPtr[idx] = oldMapPtr[i];
									break;
								}
								//move on
								idx++;
								if(idx == (newCap * 2)){
									idx = 0;
								}
							}
						}
					}
					//free old map ptr
					globalAllocator.freePow2(oldMapPtr, oldCap * 2 * sizeof(Neigh));
				}



			}







			capacity = getNextPow2(capacity * 2);
			Neigh* __restrict newPtr = (Neigh*)globalAllocator.allocPow2(capacity * sizeof(Neigh));

			if(degree <= TH0){	//Going from Type 1 to Type 2 or 3
				memcpy(newPtr, etype.type1.neigh, degree * sizeof(Neigh));
				etype.type2_3.mapArr = nullptr;
			}
			else{				//Type 2 or 3
				memcpy(newPtr, etype.type2_3.neighArr, degree * sizeof(Neigh));
				globalAllocator.freePow2(etype.type2_3.neighArr, capacity / 2 * sizeof(Neigh));
			}
			etype.type2_3.neighArr = newPtr;

			//Grow hash table if needed
			rebuildHashTable(capacity / 2, capacity);
		}

		Neigh* __restrict currNeighArr;

		if(capacity <= TH0){
			currNeighArr = etype.type1.neigh;
		}
		else{
			currNeighArr = etype.type2_3.neighArr;
		}

		//search and insert if not found
		if(capacity <= TH1){
			//Type 1 or 2, do linear search
			for(u64 i = 0; i < degree; i++){
				if(currNeighArr[i].node == dstId){
					//found same edge, just update
					currNeighArr[i].setWeight(weight);
					return;
				}
			}
		}
		else{
			//type 3, use hash table + adj list
			u32 idx = dstId & (capacity * 2 - 1);
			DstLocPair* __restrict locMap = etype.type2_3.mapArr;
			DstLocPair* __restrict insLoc = nullptr;
			//probe = 0;
			while(true){
				//probe++;
				if(locMap[idx].dst == FLAG_EMPTY_SLOT){
					//edge not found, insert
					if(insLoc){
						locMap = insLoc;	//points to the first tomb stone found
					}
					locMap[idx].dst = dstId;
					locMap[idx].loc = degree;
					break;
				}
				else if((locMap[idx].dst == FLAG_TOMB_STONE) && (insLoc == nullptr)){
					insLoc = locMap + idx;
				}
				else if(locMap[idx].dst == dstId){
					//edge found, update weight
					currNeighArr[locMap[idx].loc].setWeight(weight);
					//probingDist[probe]++;
					return;
				}
				//move on
				idx++;
				if(idx == (capacity * 2)){
					idx = 0;
				}
			}
		}
		//not found, insert
		currNeighArr[degree].node = dstId;
		currNeighArr[degree].setWeight(weight);
		degree++;
		edgeCnt++;
	}


	void deleteEdge(const Idx dstId, u64& edgeCnt){
		Neigh* __restrict currNeighArr;
		Neigh* __restrict nn = nullptr;

		if(capacity <= TH0){
			currNeighArr = etype.type1.neigh;
		}
		else{
			currNeighArr = etype.type2_3.neighArr;
		}

		//search
		if(capacity <= TH1){
			//Type 1 or 2, do linear search
			for(u64 i = 0; i < degree; i++){
				if(currNeighArr[i].node == dstId){
					nn = currNeighArr + i;
					break;
				}
			}
			if(__builtin_expect(nn != nullptr, true)){
				//edge found, delete
				degree--;
				edgeCnt--;
				nn->node = currNeighArr[degree].node;
				nn->setWeight(currNeighArr[degree].getWeight());
			}
			else{
				//edge not found, nothing to do
				return;
			}
		}
		else{
			//using hashed mode
			u32 idx = dstId & (capacity * 2 - 1);
			DstLocPair* __restrict locMap = etype.type2_3.mapArr;
			while(true){
				if(locMap[idx].dst == dstId){
					//edge found, delete
					degree--;
					edgeCnt--;
					//delSucc++;
					locMap[idx].dst = FLAG_TOMB_STONE; 				//invalidate previous hash-table entry

					const u32 loc = locMap[idx].loc;
					if(__builtin_expect(loc != degree, true)){		//nothing to do if last entry is removed
						const u32 node = currNeighArr[degree].node;
						//copy last entry
						currNeighArr[loc] = currNeighArr[degree];

						//point to correct location of the swapped entry
						u32 idxMoved = node & (capacity * 2 - 1);
						while(locMap[idxMoved].dst != node){
							idxMoved++;
							if(idxMoved == (capacity * 2)){
								idxMoved = 0;
							}
						}
						locMap[idxMoved].loc = loc;
					}
					break;
				}
				else if (locMap[idx].dst == FLAG_EMPTY_SLOT) {
					//edge not found, return
					return;
				}
				//move on
				idx++;
				if(idx == (capacity * 2)){
					idx = 0;
				}
			}
		}

		if((capacity > TH0) && ((degree * 4) <= capacity)){
			//time to reduce capacity
			const u64 oldCap = capacity;
			const u64 newCap = capacity / 2;
			capacity = newCap;

			Neigh* __restrict oldPtr = etype.type2_3.neighArr;
			Neigh* __restrict newPtr;

			if(newCap <= TH0){
				//moving from type 2 or 3 to type 1
				newPtr = etype.type1.neigh;
				capacity = TH0;
			}
			else{
				etype.type2_3.neighArr = (Neigh*)globalAllocator.allocPow2(newCap * sizeof(Neigh));
				newPtr = etype.type2_3.neighArr;
			}

			//copy old adjList and free
			memcpy(newPtr, oldPtr, degree * sizeof(Neigh));
			globalAllocator.freePow2(oldPtr, oldCap * sizeof(Neigh));

			//shrink or delete hash table if needed
			rebuildHashTable(oldCap, newCap);
		}
	}
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

#include "GraphTangoHash.h"

template <typename Neigh>
class Vertex{
public:
	 GraphTangoHash<Neigh>		inEdges;
	 GraphTangoHash<Neigh>		outEdges;
};

#endif


#ifdef USE_CAHCE_FRIENDLY_HASH_ONLY

#include "GraphTangoHash.h"

template <typename Neigh>
class Vertex{
public:
	 GraphTangoHash<Neigh>		inEdges;
	 GraphTangoHash<Neigh>		outEdges;
};

#endif


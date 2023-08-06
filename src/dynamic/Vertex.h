#pragma once

#include <cstdlib>
#include <vector>
#include "common.h"
#include "CustomAllocator.h"

#define 	FLAG_EMPTY_SLOT			0xFFFFFFFFU
#define 	FLAG_TOMB_STONE			0xFFFFFFFEU
#define		CACHE_LINE_SIZE			64

template <typename Neigh>
class alignas(CACHE_LINE_SIZE) EdgeArray{

private:
	void rebuildHashTable(u64 newCap){
		//allocate new map
		etype.type2_3.mapArr = (DstLocPair*)globalAllocator.allocate(newCap * sizeof(DstLocPair));
		etype.type2_3.mapSize = newCap;

		DstLocPair* __restrict locMap = etype.type2_3.mapArr;
		memset(locMap, -1, newCap * sizeof(DstLocPair));

		const u32 mask = newCap - 1;
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
				if(idx == newCap){
					idx = 0;
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
			u64 mapSize = 0;
		} type2_3;
	} etype;

	u8 __pad[CACHE_LINE_SIZE - sizeof(u32) - sizeof(u32) - sizeof(etype)];

	void insertEdge(const Idx dstId, const Weight weight, u64& edgeCnt){
		//First, check if needs expanding
		if(__builtin_expect(degree == capacity, false)){
			capacity = getNextPow2(capacity * 2);
			Neigh* __restrict newPtr = (Neigh*)globalAllocator.allocPow2(capacity * sizeof(Neigh));

			//handle old type
			if(degree <= TH0){
				memcpy(newPtr, etype.type1.neigh, degree * sizeof(Neigh));
			}
			else {
				memcpy(newPtr, etype.type2_3.neighArr, degree * sizeof(Neigh));
				globalAllocator.freePow2(etype.type2_3.neighArr, degree * sizeof(Neigh));
			}
			etype.type2_3.neighArr = newPtr;

			//handle new type
			if(degree <= TH1 && capacity > TH1){
				rebuildHashTable(capacity);
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
			//Grow hash table if needed
			float loadFactor = 1.0 * degree / etype.type2_3.mapSize;
			if(loadFactor >= GT_MAX_LOAD_FACTOR) {
				//free old map
				globalAllocator.freePow2(etype.type2_3.mapArr, etype.type2_3.mapSize * sizeof(DstLocPair));
				rebuildHashTable(etype.type2_3.mapSize * 2);
			}
			u32 idx = dstId & (etype.type2_3.mapSize - 1);
			DstLocPair* __restrict locMap = etype.type2_3.mapArr;
			i64 insLoc = -1LL;
			//probe = 0;
			while(true){
				//probe++;
				if(locMap[idx].dst == FLAG_EMPTY_SLOT){
					//edge not found, insert
					if(insLoc != -1LL){
						idx = insLoc;	//points to the first tomb stone found, if any
					}
					locMap[idx].dst = dstId;
					locMap[idx].loc = degree;
					break;
				}
				else if((locMap[idx].dst == FLAG_TOMB_STONE) && (insLoc == -1LL)){
					insLoc = idx;	//found a deleted slot
				}
				else if(locMap[idx].dst == dstId){
					//edge found, update weight
					currNeighArr[locMap[idx].loc].setWeight(weight);
					//probingDist[probe]++;
					return;
				}
				//move on
				idx++;
				if(idx == etype.type2_3.mapSize){
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
			u32 idx = dstId & (etype.type2_3.mapSize - 1);
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
						u32 idxMoved = node & (etype.type2_3.mapSize - 1);
						while(locMap[idxMoved].dst != node){
							idxMoved++;
							if(idxMoved == etype.type2_3.mapSize){
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
				if(idx == etype.type2_3.mapSize){
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
			//rebuildHashTable(oldCap, newCap);
		}   
	}
};


template <typename Neigh>
class Vertex{
public:
	 EdgeArray<Neigh>		inEdges;
	 EdgeArray<Neigh>		outEdges;
};

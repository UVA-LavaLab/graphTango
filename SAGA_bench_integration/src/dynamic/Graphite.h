#pragma once

#include <cstdint>
#include <stack>
#include <vector>
#include <variant>
#include <cstring>
#include <algorithm>
#include <omp.h>
#include <string>
#include <iostream>
#include <atomic>
#include <unordered_map>
#include <queue>
#include <cassert>
#include <immintrin.h>
#include <fstream>
#include <map>
#include <bitset>

#include "omp.h"

#include "abstract_data_struc.h"
#include "types.h"

#include "LockFreePoolWithList.h"
#include "CustomAllocator.h"
#include "VertexArray.h"

#include "GEmpty.h"
#include "common.h"
#include "global.h"

using namespace std;

template<typename Neigh>
class Graphite : public dataStruc {

public:

	void print(void) override {
//		std::cout << "Inserts--------------------" << std::endl;
//		std::cout << "    Total: " << insTot << std::endl;
//		std::cout << "    Succ : " << insSucc << std::endl;
//		std::cout << "    Fail : " << insTot - insSucc << std::endl;
//		std::cout << std::endl;
//
//		std::cout << "Deletes--------------------" << std::endl;
//		std::cout << "    Total: " << delTot << std::endl;
//		std::cout << "    Succ : " << delSucc << std::endl;
//		std::cout << "    Fail : " << delTot - delSucc << std::endl;
//		std::cout << std::endl;
//
//		std::cout << "Final number of edges: " << insSucc - delSucc << std::endl;
		ofstream out("probing_dist.csv");
		for(auto it : probingDist){
			out << it.first << "," << it.second << endl;
		}
	}

#ifdef USE_HYBRID_HASHMAP

	VertexArray<Neigh> vArray;
	const int num_threads;

	typedef struct{
		u64 edgeCnt = 0;
		u64 nodeCnt = 0;
		u8 pad[48];
	} ThreadInfo;

	alignas(64) ThreadInfo thInfo[32];

	Graphite(bool weighted, bool directed, i64 numNodes, i64 numThreads) : dataStruc(weighted, directed), num_threads(numThreads){
#ifdef _OPENMP
		if(numThreads > 0){
			omp_set_num_threads(numThreads);
		}
#endif

		vArray.resize(numNodes);

		#pragma omp parallel for
		for(u64 i = 0; i < numNodes; i++){
			vArray.outDegree[i] = 0;
			vArray.outCapacity[i] = 0;
			vArray.outNeighArr[i] = nullptr;
			vArray.outMapArr[i] = nullptr;

			vArray.inDegree[i] = 0;
			vArray.inCapacity[i] = 0;
			vArray.inNeighArr[i] = nullptr;
			vArray.inMapArr[i] = nullptr;
		}

		//vArray.usingHash.resize(numNodes, false);
		//vArray.dstLocMap.resize(numNodes);

		cout << "Sizeof ThreadInfo: " << sizeof(ThreadInfo) << endl;

		property.resize(numNodes, -1);
		affected.resize(numNodes);
		affected.fill(false);
	}

	~Graphite(){

	}

	int64_t in_degree(NodeID n) override {
		return vArray.inDegree[n];
	}

	int64_t out_degree(NodeID n) override {
		return vArray.outDegree[n];
	}


	void insertEdge(u64& deg, u64& cap, Neigh* __restrict &neighs, graphite_hashmap* __restrict &locMap, const Idx dstId, const Weight weight, u64& edgeCnt){
		//insTot++;
		if(deg == cap){
			const u64 newCap = getNextPow2MinRet(cap * 2);

			Neigh* __restrict oldPtr = neighs;
			Neigh* __restrict newPtr = (Neigh*)globalAllocator.allocPow2(newCap * sizeof(Neigh));

			cap = newCap;
			neighs = newPtr;

			if(deg > 0){
				memcpy(newPtr, oldPtr, deg * sizeof(Neigh));
				globalAllocator.freePow2(oldPtr, newCap / 2 * sizeof(Neigh));
			}

			if(__builtin_expect(newCap == HYBRID_HASH_PARTITION, 0)){
				//switch from linear to hashed mode
				locMap = (graphite_hashmap*)globalAllocator.allocate(sizeof(graphite_hashmap));
				new (locMap) graphite_hashmap();

				//add existing nodes to hash
				const Neigh* __restrict neighs = newPtr;
				for(u64 i = 0; i < deg; i++){
					(*locMap)[neighs[i].node] = i;
				}
			}
		}

		//search for existing edge
		if(!locMap){
			//using linear mode
			Neigh* __restrict nn = nullptr;
			for(u64 i = 0; i < deg; i++){
				if(neighs[i].node == dstId){
					nn = neighs + i;
					break;
				}
			}
			if(!nn){
				//edge not found, insert
				neighs[deg].node = dstId;
				neighs[deg].setWeight(weight);
				deg++;

				edgeCnt++;
				//insSucc++;
			}
			else{
				//edge found, update
				nn->setWeight(weight);
			}
		}
		else {
			//using hashed mode
			const auto& it = locMap->find(dstId);
			if(it == locMap->end()){
				//edge not found, insert
				Neigh& nn = neighs[deg];
				nn.node = dstId;
				nn.setWeight(weight);
				(*locMap)[dstId] = deg;
				deg++;

				edgeCnt++;
				//insSucc++;
			}
			else{
				//edge found, update weight
				neighs[it->second].setWeight(weight);
			}
		}
	}


	void deleteEdge(u64& deg, u64& cap, Neigh* __restrict &neighs, graphite_hashmap* __restrict &locMap, const Idx dstId, u64& edgeCnt){
		//delTot++;
		//search for existing edge
		if(!locMap){
			//using linear mode
			Neigh* __restrict nn = nullptr;
			for(u64 i = 0; i < deg; i++){
				if(neighs[i].node == dstId){
					nn = neighs + i;
					break;
				}
			}
			if(!nn){
				//edge not found, nothing to do
				return;
			}
			else{
				//edge found, delete
				deg--;
				edgeCnt--;
				//delSucc++;

				//move last elem
				*nn = neighs[deg];
			}
		}
		else {
			//using hashed mode
			const auto& it = locMap->find(dstId);
			if(it == locMap->end()){
				//edge not found, nothing to do
				return;
			}
			else{
				//edge found, delete
				deg--;
				edgeCnt--;
				//delSucc++;

				const u64 loc = it->second;
				locMap->erase(it);

				if(__builtin_expect(loc != deg, false)){
					//not the last element... move
					neighs[loc] = neighs[deg];
					(*locMap)[neighs[loc].node] = loc;
				}
			}
		}

		if((cap > 4) && ((deg * 4) <= cap)){
			const u64 newCap = cap / 2;

			Neigh* __restrict oldPtr = neighs;
			Neigh* __restrict newPtr = (Neigh*)globalAllocator.allocPow2(newCap * sizeof(Neigh));

			cap = newCap;
			neighs = newPtr;

			memcpy(newPtr, oldPtr, deg * sizeof(Neigh));
			globalAllocator.freePow2(oldPtr, newCap * 2 * sizeof(Neigh));

			if(__builtin_expect(newCap < HYBRID_HASH_PARTITION, 0)){
				//switch from hashed to linear mode
				globalAllocator.freePow2(locMap, sizeof(graphite_hashmap));
				locMap = nullptr;
			}
		}
	}


	void update(const EdgeList& el) override {
		const u64 batchSize = el.size();

#ifdef _OPENMP
		int thMask = (1 << getNextPow2Log2(num_threads)) - 1;

		#pragma omp parallel
		for(u64 i = 0; i < batchSize; i++){
			const i64 src = el[i].source;
			const i64 dst = el[i].destination;
			const i64 actualTh = omp_get_thread_num();

			i64 targetTh = (src / 64) & thMask;
			if(targetTh == actualTh){
				if(!el[i].sourceExists){
					thInfo[actualTh].nodeCnt++;
				}
				if(!el[i].destExists){
					thInfo[actualTh].nodeCnt++;
				}

				if(!affected[src]){
					affected[src] = true;
				}

				if(!el[i].isDelete){
					//insert out edge
					insertEdge(vArray.outDegree[src], vArray.outCapacity[src], vArray.outNeighArr[src], vArray.outMapArr[src], dst, el[i].weight, thInfo[actualTh].edgeCnt);
				}
				else{
					//delete out edge
					deleteEdge(vArray.outDegree[src], vArray.outCapacity[src], vArray.outNeighArr[src], vArray.outMapArr[src], dst, thInfo[actualTh].edgeCnt);
				}
			}

			targetTh = (dst / 64) & thMask;
			if(targetTh == actualTh){
				if(!affected[dst]){
					affected[dst] = true;
				}

				u64 garbage;
				if(!el[i].isDelete){
					//insert in edge
					insertEdge(vArray.inDegree[dst], vArray.inCapacity[dst], vArray.inNeighArr[dst], vArray.inMapArr[dst], src, el[i].weight, garbage);
				}
				else{
					//delete in edge
					deleteEdge(vArray.inDegree[dst], vArray.inCapacity[dst], vArray.inNeighArr[dst], vArray.inMapArr[dst], src, garbage);
				}
			}

		}
#else
		for(u64 i = 0; i < batchSize; i++){
			const u64 src = el[i].source;
			const u64 dst = el[i].destination;

			if(!el[i].sourceExists){
				thInfo[0].nodeCnt++;
			}
			if(!el[i].destExists){
				thInfo[0].nodeCnt++;
			}
			//we do not need atomic operation on affected as long as "some" thread updates it
			if(!affected[src]){
				affected[src] = true;
			}
			if(!affected[dst]){
				affected[dst] = true;
			}

			u64 garbage;
			if(!el[i].isDelete){
				//insert out edge
				insertEdge(vArray.outDegree[src], vArray.outCapacity[src], vArray.outNeighArr[src], vArray.outMapArr[src], dst, el[i].weight, thInfo[0].edgeCnt);

				//insert in edge
				insertEdge(vArray.inDegree[dst], vArray.inCapacity[dst], vArray.inNeighArr[dst], vArray.inMapArr[dst], src, el[i].weight, garbage);
			}
			else{
				//delete out edge
				deleteEdge(vArray.outDegree[src], vArray.outCapacity[src], vArray.outNeighArr[src], vArray.outMapArr[src], dst, thInfo[0].edgeCnt);

				//delete in edge
				deleteEdge(vArray.inDegree[dst], vArray.inCapacity[dst], vArray.inNeighArr[dst], vArray.inMapArr[dst], src, garbage);
			}
		}
#endif

		for(u64 i = 0; i < num_threads; i++){
			num_edges += thInfo[i].edgeCnt;
			thInfo[i].edgeCnt = 0;
			num_nodes += thInfo[i].nodeCnt;
			thInfo[i].nodeCnt = 0;
		}
		//num_nodes = el[batchSize - 1].lastAssignedId + 1;
	}

#endif



#ifdef USE_HYBRID_HASHMAP_WITH_CFH

#define FLAG_EMPTY_SLOT			0xFFFFFFFFU
#define FLAG_TOMB_STONE			0xFFFFFFFEU

	VertexArray<Neigh> vArray;
	const int num_threads;

	map<u64, u64> probingDist;
	u64 probe;

	typedef struct{
		u64 edgeCnt = 0;
		u64 nodeCnt = 0;
		u8 pad[48];
	} ThreadInfo;

	alignas(64) ThreadInfo thInfo[32];

	Graphite(bool weighted, bool directed, i64 numNodes, i64 numThreads) : dataStruc(weighted, directed), num_threads(numThreads){
#ifdef _OPENMP
		if(numThreads > 0){
			omp_set_num_threads(numThreads);
		}
#endif

		vArray.resize(numNodes);

		#pragma omp parallel for
		for(u64 i = 0; i < numNodes; i++){
			vArray.outDegree[i] = 0;
			vArray.outCapacity[i] = 0;
			vArray.outNeighArr[i] = nullptr;
			vArray.outMapArr[i] = nullptr;

			vArray.inDegree[i] = 0;
			vArray.inCapacity[i] = 0;
			vArray.inNeighArr[i] = nullptr;
			vArray.inMapArr[i] = nullptr;
		}

		//vArray.usingHash.resize(numNodes, false);
		//vArray.dstLocMap.resize(numNodes);

		cout << "Sizeof ThreadInfo: " << sizeof(ThreadInfo) << endl;

		property.resize(numNodes, -1);
		affected.resize(numNodes);
		affected.fill(false);
	}

	~Graphite(){

	}

	int64_t in_degree(NodeID n) override {
		return vArray.inDegree[n];
	}

	int64_t out_degree(NodeID n) override {
		return vArray.outDegree[n];
	}

	//for now, use linear probing
	inline DstLocPair* findInsertionPoint(DstLocPair* __restrict &locMap, u32 dst, u32 tableSize){
		u32 idx = dst & (tableSize - 1);
		while(true){
			if((locMap[idx].dst == dst) || (locMap[idx].dst == FLAG_EMPTY_SLOT)){
				//found key or an empty slot
				return locMap + idx;
			}
			//move on
			idx++;
			if(idx == tableSize){
				idx = 0;
			}
		}
		return nullptr;
	}

	void insertEdge(u64& deg, u64& cap, Neigh* __restrict &neighs, DstLocPair* __restrict &locMap, const Idx dstId, const Weight weight, u64& edgeCnt){
		//insTot++;
		if(deg == cap){
			const u64 newCap = getNextPow2MinRet(cap * 2);

			Neigh* __restrict oldPtr = neighs;
			Neigh* __restrict newPtr = (Neigh*)globalAllocator.allocPow2(newCap * sizeof(Neigh));

			cap = newCap;
			neighs = newPtr;

			if(deg > 0){
				memcpy(newPtr, oldPtr, deg * sizeof(Neigh));
				globalAllocator.freePow2(oldPtr, newCap / 2 * sizeof(Neigh));
			}

			if(newCap >= HYBRID_HASH_PARTITION){
				if(locMap){
					//free old loc map
					globalAllocator.freePow2(locMap, cap * sizeof(DstLocPair));
				}
				locMap = (DstLocPair*)globalAllocator.allocate(sizeof(DstLocPair) * cap * 2);
				memset(locMap, -1, sizeof(DstLocPair) * cap * 2);

				const u32 mask = cap * 2 - 1;

				//add existing nodes to hash
				const Neigh* __restrict nn = neighs;
				for(u64 i = 0; i < deg; i++){
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
						if(idx == (cap * 2)){
							idx = 0;
						}
					}
				}
			}
		}

		//search for existing edge
		if(!locMap){
			//using linear mode
			Neigh* __restrict nn = nullptr;
			for(u64 i = 0; i < deg; i++){
				if(neighs[i].node == dstId){
					nn = neighs + i;
					break;
				}
			}
			if(!nn){
				//edge not found, insert
				neighs[deg].node = dstId;
				neighs[deg].setWeight(weight);
				deg++;

				edgeCnt++;
				//insSucc++;
			}
			else{
				//edge found, update
				nn->setWeight(weight);
			}
		}
		else {
			//using hashed mode
			u32 idx = dstId & (cap * 2 - 1);
			probe = 0;
			while(true){
				probe++;
				if(locMap[idx].dst == FLAG_EMPTY_SLOT){

					//edge not found, insert
					Neigh& nn = neighs[deg];
					nn.node = dstId;
					nn.setWeight(weight);
					locMap[idx].dst = dstId;
					locMap[idx].loc = deg;
					deg++;

					edgeCnt++;
					//insSucc++;
					probingDist[probe]++;
					return;
				}
				else if(locMap[idx].dst == dstId){
					//edge found, update weight
					neighs[locMap[idx].loc].setWeight(weight);
					probingDist[probe]++;
					return;
				}
				//move on
				idx++;
				if(idx == (cap * 2)){
					idx = 0;
				}
			}
		}
	}



	void deleteEdge(u64& deg, u64& cap, Neigh* __restrict &neighs, DstLocPair* __restrict &locMap, const Idx dstId, u64& edgeCnt){
		//delTot++;
		//search for existing edge
		if(!locMap){
			//using linear mode
			Neigh* __restrict nn = nullptr;
			for(u64 i = 0; i < deg; i++){
				if(neighs[i].node == dstId){
					nn = neighs + i;
					break;
				}
			}
			if(!nn){
				//edge not found, return
				return;
			}
			else{
				//edge found, delete
				deg--;
				edgeCnt--;
				//delSucc++;
				nn->node = neighs[deg].node;
				nn->setWeight(neighs[deg].getWeight());
			}
		}
		else {
			//using hashed mode
			u32 idx = dstId & (cap * 2 - 1);
			while(true){
				if(locMap[idx].dst == FLAG_EMPTY_SLOT){
					//edge not found, return
					return;
				}
				else if(locMap[idx].dst == dstId){
					//edge found, delete
					deg--;
					edgeCnt--;
					//delSucc++;
					locMap[idx].dst = FLAG_TOMB_STONE; 				//invalidate previous hash-table entry

					const u32 loc = locMap[idx].loc;
					if(__builtin_expect(loc != deg, true)){		//nothing to do if last entry is removed
						const u32 node = neighs[deg].node;
						//copy last entry
						neighs[loc] = neighs[deg];

						//point to correct location of the swapped entry
						u32 idxMoved = node & (cap * 2 - 1);
						while(locMap[idxMoved].dst != node){
							idxMoved++;
							if(idxMoved == (cap * 2)){
								idxMoved = 0;
							}
						}
						locMap[idxMoved].loc = loc;
					}
					break;
				}
				//move on
				idx++;
				if(idx == (cap * 2)){
					idx = 0;
				}
			}
		}

		if((cap > 4) && (deg * 4) <= cap){
			//time to reduce capacity

			const u64 newCap = cap / 2;
			Neigh* __restrict oldPtr = neighs;
			Neigh* __restrict newPtr = (Neigh*)globalAllocator.allocPow2(newCap * sizeof(Neigh));

			//copy old adjList and free
			memcpy(newPtr, oldPtr, deg * sizeof(Neigh));
			globalAllocator.freePow2(oldPtr, cap * sizeof(Neigh));

			cap = newCap;
			neighs = newPtr;

			if(locMap){
				//free old loc map (if any)
				globalAllocator.freePow2(locMap, cap * sizeof(DstLocPair));
				locMap = nullptr;
			}

			if(cap >= HYBRID_HASH_PARTITION){
				locMap = (DstLocPair*)globalAllocator.allocate(sizeof(DstLocPair) * cap * 2);
				memset(locMap, -1, sizeof(DstLocPair) * cap * 2);	//This is bad. Is there any way to circumvent this?

				const u32 mask = cap * 2 - 1;

				//add existing nodes to hash
				const Neigh* __restrict nn = neighs;
				for(u64 i = 0; i < deg; i++){
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
						if(idx == (cap * 2)){
							idx = 0;
						}
					}
				}
			}
		}
	}

	void update(const EdgeList& el) override {
		probe = 0;
		const u64 batchSize = el.size();

#ifdef _OPENMP
		int thMask = (1 << getNextPow2Log2(num_threads)) - 1;

		#pragma omp parallel
		for(u64 i = 0; i < batchSize; i++){
			const i64 src = el[i].source;
			const i64 dst = el[i].destination;
			const i64 actualTh = omp_get_thread_num();

			i64 targetTh = (src / 64) & thMask;
			if(targetTh == actualTh){
				if(!el[i].sourceExists){
					thInfo[actualTh].nodeCnt++;
				}
				if(!el[i].destExists){
					thInfo[actualTh].nodeCnt++;
				}

				if(!affected[src]){
					affected[src] = true;
				}

				if(!el[i].isDelete){
					//insert out edge
					insertEdge(vArray.outDegree[src], vArray.outCapacity[src], vArray.outNeighArr[src], vArray.outMapArr[src], dst, el[i].weight, thInfo[actualTh].edgeCnt);
				}
				else{
					//delete out edge
					deleteEdge(vArray.outDegree[src], vArray.outCapacity[src], vArray.outNeighArr[src], vArray.outMapArr[src], dst, thInfo[0].edgeCnt);
				}
			}

			targetTh = (dst / 64) & thMask;
			if(targetTh == actualTh){
				if(!affected[dst]){
					affected[dst] = true;
				}

				u64 garbage;
				if(!el[i].isDelete){
					//insert in edge
					insertEdge(vArray.inDegree[dst], vArray.inCapacity[dst], vArray.inNeighArr[dst], vArray.inMapArr[dst], src, el[i].weight, garbage);
				}
				else{
					//delete in edge
					deleteEdge(vArray.inDegree[dst], vArray.inCapacity[dst], vArray.inNeighArr[dst], vArray.inMapArr[dst], src, garbage);
				}
			}

		}
#else
		for(u64 i = 0; i < batchSize; i++){
			const u64 src = el[i].source;
			const u64 dst = el[i].destination;

			if(!el[i].sourceExists){
				thInfo[0].nodeCnt++;
			}
			if(!el[i].destExists){
				thInfo[0].nodeCnt++;
			}
			//we do not need atomic operation on affected as long as "some" thread updates it
			if(!affected[src]){
				affected[src] = true;
			}
			if(!affected[dst]){
				affected[dst] = true;
			}

			u64 garbage;
			if(!el[i].isDelete){
				//insertion
				//insert out edge
				insertEdge(vArray.outDegree[src], vArray.outCapacity[src], vArray.outNeighArr[src], vArray.outMapArr[src], dst, el[i].weight, thInfo[0].edgeCnt);

				//insert in edge
				insertEdge(vArray.inDegree[dst], vArray.inCapacity[dst], vArray.inNeighArr[dst], vArray.inMapArr[dst], src, el[i].weight, garbage);
			}
			else{
				//delete out edge
				deleteEdge(vArray.outDegree[src], vArray.outCapacity[src], vArray.outNeighArr[src], vArray.outMapArr[src], dst, thInfo[0].edgeCnt);

				//delete in edge
				deleteEdge(vArray.inDegree[dst], vArray.inCapacity[dst], vArray.inNeighArr[dst], vArray.inMapArr[dst], src, garbage);
			}
		}
#endif

		for(u64 i = 0; i < num_threads; i++){
			num_edges += thInfo[i].edgeCnt;
			thInfo[i].edgeCnt = 0;
			num_nodes += thInfo[i].nodeCnt;
			thInfo[i].nodeCnt = 0;
		}
		//num_nodes = el[batchSize - 1].lastAssignedId + 1;
	}

#endif



#ifdef USE_HYBRID_HASHMAP_WITH_GROUPING

	Vertex<Neigh>* vArray;
	const int num_threads;

	typedef struct{
		u64 edgeCnt = 0;
		u64 nodeCnt = 0;
		u8 pad[48];
	} ThreadInfo;

	alignas(64) ThreadInfo thInfo[32];

	Graphite(bool weighted, bool directed, i64 numNodes, i64 numThreads) : dataStruc(weighted, directed), num_threads(numThreads){
#ifdef _OPENMP
		if(numThreads > 0){
			omp_set_num_threads(numThreads);
		}
#endif
		vArray = (Vertex<Neigh>*)aligned_alloc(64, numNodes * sizeof(Vertex<Neigh>));
		memset(vArray, 0, numNodes * sizeof(Vertex<Neigh>));

		cout << "Size of Vertex: " << sizeof(Vertex<Neigh>) << endl;
		cout << "Size of ThreadInfo: " << sizeof(ThreadInfo) << endl;

		property.resize(numNodes, -1);
		affected.resize(numNodes);
		affected.fill(false);
	}

	~Graphite(){
		free(vArray);
	}

	void print() override {
		cout << "Done" << endl;
	}

	int64_t in_degree(NodeID n) override {
		return vArray[n].inEdges.degree;
	}

	int64_t out_degree(NodeID n) override {
		return vArray[n].outEdges.degree;
	}


	void insertEdge(EdgeArray<Neigh>* __restrict edgePtr, const Idx dstId, const Weight weight, u64& edgeCnt){
		if(edgePtr->degree == edgePtr->capacity){
			const u64 newCap = getNextPow2MinRet(edgePtr->capacity * 2);

			Neigh* __restrict oldPtr = edgePtr->neighArr;
			Neigh* __restrict newPtr = (Neigh*)globalAllocator.allocPow2(newCap * sizeof(Neigh));

			edgePtr->capacity = newCap;
			edgePtr->neighArr = newPtr;

			if(edgePtr->degree > 0){
				memcpy(newPtr, oldPtr, edgePtr->degree * sizeof(Neigh));
				globalAllocator.freePow2(oldPtr, newCap / 2 * sizeof(Neigh));
			}

			if(__builtin_expect(newCap == HYBRID_HASH_PARTITION, 0)){
				//switch from linear to hashed mode
				edgePtr->mapArr = (graphite_hashmap*)globalAllocator.allocate(sizeof(graphite_hashmap));
				new (edgePtr->mapArr) graphite_hashmap();

				//add existing nodes to hash
				const Neigh* __restrict neighs = newPtr;
				for(u64 i = 0; i < edgePtr->degree; i++){
					(*edgePtr->mapArr)[neighs[i].node] = i;
				}
			}
		}

		//search for existing edge
		if(!edgePtr->mapArr){
			//using linear mode
			Neigh* __restrict nn = nullptr;
			for(u64 i = 0; i < edgePtr->degree; i++){
				if(edgePtr->neighArr[i].node == dstId){
					nn = edgePtr->neighArr + i;
					break;
				}
			}
			if(!nn){
				//edge not found, insert
				edgePtr->neighArr[edgePtr->degree].node = dstId;
				edgePtr->neighArr[edgePtr->degree].setWeight(weight);
				edgePtr->degree++;

				edgeCnt++;
			}
			else{
				//edge found, update
				nn->setWeight(weight);
			}
		}
		else {
			//using hashed mode
			const auto& it = edgePtr->mapArr->find(dstId);
			if(it == edgePtr->mapArr->end()){
				//edge not found, insert
				Neigh& nn = edgePtr->neighArr[edgePtr->degree];
				nn.node = dstId;
				nn.setWeight(weight);
				(*edgePtr->mapArr)[dstId] = edgePtr->degree;
				edgePtr->degree++;

				edgeCnt++;
			}
			else{
				//edge found, update weight
				edgePtr->neighArr[it->second].setWeight(weight);
			}
		}
	}



	void update(const EdgeList& el) override {
		const u64 batchSize = el.size();

#ifdef _OPENMP
		int thMask = (1 << getNextPow2Log2(num_threads)) - 1;

		#pragma omp parallel
		for(u64 i = 0; i < batchSize; i++){
			const i64 src = el[i].source;
			const i64 dst = el[i].destination;
			const i64 actualTh = omp_get_thread_num();

			i64 targetTh = (src / 64) & thMask;
			if(targetTh == actualTh){
				if(!el[i].sourceExists){
					thInfo[actualTh].nodeCnt++;
				}
				if(!el[i].destExists){
					thInfo[actualTh].nodeCnt++;
				}

				if(!affected[src]){
					affected[src] = true;
				}

				//insert out edge
				insertEdge(&vArray[src].outEdges, dst, el[i].weight, thInfo[actualTh].edgeCnt);
			}

			targetTh = (dst / 64) & thMask;
			if(targetTh == actualTh){
				if(!affected[dst]){
					affected[dst] = true;
				}
				//insert in edge
				insertEdge(&vArray[dst].inEdges, src, el[i].weight, thInfo[actualTh].edgeCnt);
			}

		}
#else
		for(u64 i = 0; i < batchSize; i++){
			const u64 src = el[i].source;
			const u64 dst = el[i].destination;

			if(!el[i].sourceExists){
				thInfo[0].nodeCnt++;
			}
			if(!el[i].destExists){
				thInfo[0].nodeCnt++;
			}
			if(!affected[src]){
				affected[src] = true;
			}
			if(!affected[dst]){
				affected[dst] = true;
			}

			//insert out edge
			insertEdge(&vArray[src].outEdges, dst, el[i].weight, thInfo[0].edgeCnt);

			//insert in edge
			u64 garbage;
			insertEdge(&vArray[dst].inEdges, src, el[i].weight, garbage);
		}
#endif

		for(u64 i = 0; i < num_threads; i++){
			num_edges += thInfo[i].edgeCnt;
			thInfo[i].edgeCnt = 0;
			num_nodes += thInfo[i].nodeCnt;
			thInfo[i].nodeCnt = 0;
		}
		//num_nodes = el[batchSize - 1].lastAssignedId + 1;
	}
#endif


#ifdef USE_SORTED_EDGES

	Vertex<Neigh>* vArray;
	const int num_threads;

	typedef struct{
		u64 edgeCnt = 0;
		u64 nodeCnt = 0;
		u8 pad[48];
	} ThreadInfo;

	alignas(64) ThreadInfo thInfo[32];

	Graphite(bool weighted, bool directed, i64 numNodes, i64 numThreads) : dataStruc(weighted, directed), num_threads(numThreads){
#ifdef _OPENMP
		if(numThreads > 0){
			omp_set_num_threads(numThreads);
		}
#endif
		vArray = (Vertex<Neigh>*)aligned_alloc(64, numNodes * sizeof(Vertex<Neigh>));
		memset(vArray, 0, numNodes * sizeof(Vertex<Neigh>));

		cout << "Size of Vertex: " << sizeof(Vertex<Neigh>) << endl;
		cout << "Size of ThreadInfo: " << sizeof(ThreadInfo) << endl;

		property.resize(numNodes, -1);
		affected.resize(numNodes);
		affected.fill(false);
	}

	~Graphite(){
		free(vArray);
	}

	void print() override {
		cout << "Done" << endl;
	}

	int64_t in_degree(NodeID n) override {
		return vArray[n].inEdges.degree;
	}

	int64_t out_degree(NodeID n) override {
		return vArray[n].outEdges.degree;
	}


	void insertEdge(EdgeArray<Neigh>* __restrict edgePtr, const Idx dstId, const Weight weight, u64& edgeCnt){
		Neigh* __restrict insLoc = nullptr;

		//first, do linear search
		for(u64 i = 0; i < edgePtr->degree; i++){
			if(edgePtr->neighArr[i].node == dstId){
				//found
				insLoc = edgePtr->neighArr + i;
				break;
			}
		}

		if(!insLoc && edgePtr->degree > LINEAR_BUFF_SIZE){
			//do binary search
			u64 low = LINEAR_BUFF_SIZE;
			u64 high = edgePtr->degree;
		}

		if(!insLoc){
			//not found, insert
			if(edgePtr->degree == edgePtr->capacity){
				const u64 newCap = getNextPow2MinRet(edgePtr->capacity * 2);

				Neigh* __restrict oldPtr = edgePtr->neighArr;
				Neigh* __restrict newPtr = (Neigh*)globalAllocator.allocPow2(newCap * sizeof(Neigh));

				edgePtr->capacity = newCap;
				edgePtr->neighArr = newPtr;

				if(edgePtr->degree > 0){
					memcpy(newPtr, oldPtr, edgePtr->degree * sizeof(Neigh));
					globalAllocator.freePow2(oldPtr, newCap / 2 * sizeof(Neigh));
				}
			}

			insLoc->node = dstId;
			edgePtr->degree++;

			edgeCnt++;
		}
		insLoc->setWeight(weight);
	}



	void update(const EdgeList& el) override {
		const u64 batchSize = el.size();

#ifdef _OPENMP
		int thMask = (1 << getNextPow2Log2(num_threads)) - 1;

		#pragma omp parallel
		for(u64 i = 0; i < batchSize; i++){
			const i64 src = el[i].source;
			const i64 dst = el[i].destination;
			const i64 actualTh = omp_get_thread_num();

			i64 targetTh = (src / 64) & thMask;
			if(targetTh == actualTh){
				if(!el[i].sourceExists){
					thInfo[actualTh].nodeCnt++;
				}
				if(!el[i].destExists){
					thInfo[actualTh].nodeCnt++;
				}

				if(!affected[src]){
					affected[src] = true;
				}

				//insert out edge
				insertEdge(&vArray[src].outEdges, dst, el[i].weight, thInfo[actualTh].edgeCnt);
			}

			targetTh = (dst / 64) & thMask;
			if(targetTh == actualTh){
				if(!affected[dst]){
					affected[dst] = true;
				}
				//insert in edge
				u64 garbase;
				insertEdge(&vArray[dst].inEdges, src, el[i].weight, garbase);
			}

		}
#else
		for(u64 i = 0; i < batchSize; i++){
			const u64 src = el[i].source;
			const u64 dst = el[i].destination;

			if(!el[i].sourceExists){
				thInfo[0].nodeCnt++;
			}
			if(!el[i].destExists){
				thInfo[0].nodeCnt++;
			}
			if(!affected[src]){
				affected[src] = true;
			}
			if(!affected[dst]){
				affected[dst] = true;
			}

			//insert out edge
			insertEdge(&vArray[src].outEdges, dst, el[i].weight, thInfo[0].edgeCnt);

			//insert in edge
			u64 garbage;
			insertEdge(&vArray[dst].inEdges, src, el[i].weight, garbage);
		}
#endif

		for(u64 i = 0; i < num_threads; i++){
			num_edges += thInfo[i].edgeCnt;
			thInfo[i].edgeCnt = 0;
			num_nodes += thInfo[i].nodeCnt;
			thInfo[i].nodeCnt = 0;
		}
		//num_nodes = el[batchSize - 1].lastAssignedId + 1;
	}
#endif


#ifdef USE_HYBRID_HASHMAP_WITH_GROUPING_TIGHTER

	Vertex<Neigh>* vArray;
	const int num_threads;

	typedef struct{
		u64 edgeCnt = 0;
		u64 nodeCnt = 0;
		u8 pad[48];
	} ThreadInfo;

	alignas(64) ThreadInfo thInfo[32];

	Graphite(bool weighted, bool directed, i64 numNodes, i64 numThreads) : dataStruc(weighted, directed), num_threads(numThreads){
#ifdef _OPENMP
		if(numThreads > 0){
			omp_set_num_threads(numThreads);
		}
#endif
		vArray = (Vertex<Neigh>*)aligned_alloc(64, numNodes * sizeof(Vertex<Neigh>));
		for(u64 i = 0; i < numNodes; i++){
			Vertex<Neigh>& v = vArray[i];
			v.inEdges.degree = 0;
			v.inEdges.capacity = INITIAL_EDGES;
			v.inEdges.neighArr = v.inEdges.neigh;
			v.inEdges.mapArr = nullptr;

			v.outEdges.degree = 0;
			v.outEdges.capacity = INITIAL_EDGES;
			v.outEdges.neighArr = v.outEdges.neigh;
			v.outEdges.mapArr = nullptr;
		}

		cout << "Size of EdgeArray: " << sizeof(EdgeArray<Neigh>) << endl;
		cout << "Size of Vertex: " << sizeof(Vertex<Neigh>) << endl;
		cout << "Size of ThreadInfo: " << sizeof(ThreadInfo) << endl;
		cout << "Size of Neigh: " << sizeof(Neigh) << endl;

		property.resize(numNodes, -1);
		affected.resize(numNodes);
		affected.fill(false);
	}

	~Graphite(){
		free(vArray);
	}

	void print() override {
		cout << "Done" << endl;
	}

	int64_t in_degree(NodeID n) override {
		return vArray[n].inEdges.degree;
	}

	int64_t out_degree(NodeID n) override {
		return vArray[n].outEdges.degree;
	}


	void insertEdge(EdgeArray<Neigh>* __restrict edgePtr, const Idx dstId, const Weight weight, u64& edgeCnt){
		if(edgePtr->degree == edgePtr->capacity){
			const u64 newCap = getNextPow2MinRet(edgePtr->capacity * 2);

			Neigh* __restrict oldPtr = edgePtr->neighArr;
			Neigh* __restrict newPtr = (Neigh*)globalAllocator.allocPow2(newCap * sizeof(Neigh));

			edgePtr->capacity = newCap;
			edgePtr->neighArr = newPtr;

			if(edgePtr->degree > 0){
				memcpy(newPtr, oldPtr, edgePtr->degree * sizeof(Neigh));
				if(oldPtr != edgePtr->neigh){
					globalAllocator.freePow2(oldPtr, newCap / 2 * sizeof(Neigh));
				}
			}

			if(__builtin_expect(newCap == HYBRID_HASH_PARTITION, 0)){
				//switch from linear to hashed mode
				edgePtr->mapArr = (graphite_hashmap*)globalAllocator.allocate(sizeof(graphite_hashmap));
				new (edgePtr->mapArr) graphite_hashmap();

				//add existing nodes to hash
				const Neigh* __restrict neighs = newPtr;
				for(u64 i = 0; i < edgePtr->degree; i++){
					(*edgePtr->mapArr)[neighs[i].node] = i;
				}
			}
		}

		//search for existing edge
		if(!edgePtr->mapArr){
			//using linear mode
			Neigh* __restrict nn = nullptr;
			for(u64 i = 0; i < edgePtr->degree; i++){
				if(edgePtr->neighArr[i].node == dstId){
					nn = edgePtr->neighArr + i;
					break;
				}
			}
			if(!nn){
				//edge not found, insert
				edgePtr->neighArr[edgePtr->degree].node = dstId;
				edgePtr->neighArr[edgePtr->degree].setWeight(weight);
				edgePtr->degree++;

				edgeCnt++;
			}
			else{
				//edge found, update
				nn->setWeight(weight);
			}
		}
		else {
			//using hashed mode
			const auto& it = edgePtr->mapArr->find(dstId);
			if(it == edgePtr->mapArr->end()){
				//edge not found, insert
				Neigh& nn = edgePtr->neighArr[edgePtr->degree];
				nn.node = dstId;
				nn.setWeight(weight);
				(*edgePtr->mapArr)[dstId] = edgePtr->degree;
				edgePtr->degree++;

				edgeCnt++;
			}
			else{
				//edge found, update weight
				edgePtr->neighArr[it->second].setWeight(weight);
			}
		}
	}



	void update(const EdgeList& el) override {
		const u64 batchSize = el.size();

#ifdef _OPENMP
		int thMask = (1 << getNextPow2Log2(num_threads)) - 1;

		#pragma omp parallel
		for(u64 i = 0; i < batchSize; i++){
			const i64 src = el[i].source;
			const i64 dst = el[i].destination;
			const i64 actualTh = omp_get_thread_num();

			i64 targetTh = (src / 64) & thMask;
			if(targetTh == actualTh){
				if(!el[i].sourceExists){
					thInfo[actualTh].nodeCnt++;
				}
				if(!el[i].destExists){
					thInfo[actualTh].nodeCnt++;
				}

				if(!affected[src]){
					affected[src] = true;
				}

				//insert out edge
				insertEdge(&vArray[src].outEdges, dst, el[i].weight, thInfo[actualTh].edgeCnt);
			}

			targetTh = (dst / 64) & thMask;
			if(targetTh == actualTh){
				if(!affected[dst]){
					affected[dst] = true;
				}
				//insert in edge
				u64 garbase;
				insertEdge(&vArray[dst].inEdges, src, el[i].weight, garbase);
			}

		}
#else
		for(u64 i = 0; i < batchSize; i++){
			const u64 src = el[i].source;
			const u64 dst = el[i].destination;

			if(!el[i].sourceExists){
				thInfo[0].nodeCnt++;
			}
			if(!el[i].destExists){
				thInfo[0].nodeCnt++;
			}
			if(!affected[src]){
				affected[src] = true;
			}
			if(!affected[dst]){
				affected[dst] = true;
			}

			//insert out edge
			insertEdge(&vArray[src].outEdges, dst, el[i].weight, thInfo[0].edgeCnt);

			//insert in edge
			u64 garbage;
			insertEdge(&vArray[dst].inEdges, src, el[i].weight, garbage);
		}
#endif

		for(u64 i = 0; i < num_threads; i++){
			num_edges += thInfo[i].edgeCnt;
			thInfo[i].edgeCnt = 0;
			num_nodes += thInfo[i].nodeCnt;
			thInfo[i].nodeCnt = 0;
		}
		//num_nodes = el[batchSize - 1].lastAssignedId + 1;
	}
#endif


#ifdef USE_HYBRID_HASHMAP_WITH_GROUPING_AND_EDGE_ARR_LOCKING

	Vertex<Neigh>* vArray;
	const int num_threads;

	typedef struct{
		u64 edgeCnt = 0;
		u64 nodeCnt = 0;
		u8 pad[48];
	} ThreadInfo;

	alignas(64) ThreadInfo thInfo[32];

	Graphite(bool weighted, bool directed, i64 numNodes, i64 numThreads) : dataStruc(weighted, directed), num_threads(numThreads){
#ifdef _OPENMP
		if(numThreads > 0){
			omp_set_num_threads(numThreads);
		}
#endif
		vArray = (Vertex<Neigh>*)aligned_alloc(64, numNodes * sizeof(Vertex<Neigh>));
		//memset(vArray, 0, numNodes * sizeof(Vertex<Neigh>));
		for(u64 i = 0; i < numNodes; i++){
			vArray[i].inEdges.degree = 0;
			vArray[i].inEdges.capacity = 0;
			vArray[i].inEdges.neighArr = nullptr;
			vArray[i].inEdges.mapArr = nullptr;
			vArray[i].outEdges.degree = 0;
			vArray[i].outEdges.capacity = 0;
			vArray[i].outEdges.neighArr = nullptr;
			vArray[i].outEdges.mapArr = nullptr;
#ifdef _OPENMP
			omp_init_lock(&vArray[i].inEdges.lock);
			omp_init_lock(&vArray[i].outEdges.lock);
#endif
		}

		cout << "Size of Vertex: " << sizeof(Vertex<Neigh>) << endl;
		cout << "Size of ThreadInfo: " << sizeof(ThreadInfo) << endl;

		property.resize(numNodes, -1);
		affected.resize(numNodes);
		affected.fill(false);
	}

	~Graphite(){
		free(vArray);
	}

	void print() override {
		cout << "Done" << endl;
	}

	int64_t in_degree(NodeID n) override {
		return vArray[n].inEdges.degree;
	}

	int64_t out_degree(NodeID n) override {
		return vArray[n].outEdges.degree;
	}


	void insertEdge(EdgeArray<Neigh>* __restrict edgePtr, const Idx dstId, const Weight weight, u64& edgeCnt){
		if(edgePtr->degree == edgePtr->capacity){
			const u64 newCap = getNextPow2MinRet(edgePtr->capacity * 2);

			Neigh* __restrict oldPtr = edgePtr->neighArr;
			Neigh* __restrict newPtr = (Neigh*)globalAllocator.allocPow2(newCap * sizeof(Neigh));

			edgePtr->capacity = newCap;
			edgePtr->neighArr = newPtr;

			if(edgePtr->degree > 0){
				memcpy(newPtr, oldPtr, edgePtr->degree * sizeof(Neigh));
				globalAllocator.freePow2(oldPtr, newCap / 2 * sizeof(Neigh));
			}

			if(__builtin_expect(newCap == HYBRID_HASH_PARTITION, 0)){
				//switch from linear to hashed mode
				edgePtr->mapArr = (graphite_hashmap*)globalAllocator.allocate(sizeof(graphite_hashmap));
				new (edgePtr->mapArr) graphite_hashmap();

				//add existing nodes to hash
				const Neigh* __restrict neighs = newPtr;
				for(u64 i = 0; i < edgePtr->degree; i++){
					(*edgePtr->mapArr)[neighs[i].node] = i;
				}
			}
		}

		//search for existing edge
		if(!edgePtr->mapArr){
			//using linear mode
			Neigh* __restrict nn = nullptr;
			for(u64 i = 0; i < edgePtr->degree; i++){
				if(edgePtr->neighArr[i].node == dstId){
					nn = edgePtr->neighArr + i;
					break;
				}
			}
			if(!nn){
				//edge not found, insert
				edgePtr->neighArr[edgePtr->degree].node = dstId;
				edgePtr->neighArr[edgePtr->degree].setWeight(weight);
				edgePtr->degree++;

				edgeCnt++;
			}
			else{
				//edge found, update
				nn->setWeight(weight);
			}
		}
		else {
			//using hashed mode
			const auto& it = edgePtr->mapArr->find(dstId);
			if(it == edgePtr->mapArr->end()){
				//edge not found, insert
				Neigh& nn = edgePtr->neighArr[edgePtr->degree];
				nn.node = dstId;
				nn.setWeight(weight);
				(*edgePtr->mapArr)[dstId] = edgePtr->degree;
				edgePtr->degree++;

				edgeCnt++;
			}
			else{
				//edge found, update weight
				edgePtr->neighArr[it->second].setWeight(weight);
			}
		}
	}



	void update(const EdgeList& el) override {
		const u64 batchSize = el.size();

		#pragma omp parallel for
		for(u64 i = 0; i < batchSize; i++){
			const Edge e = el[i];
			const int th = omp_get_thread_num();

			if(!e.sourceExists){
				thInfo[th].nodeCnt++;
			}
			if(!e.destExists){
				thInfo[th].nodeCnt++;
			}

			if(!affected[e.source]){
				affected[e.source] = true;
			}
			if(!affected[e.destination]){
				affected[e.destination] = true;
			}

			//take a lock on the out edge of src
			omp_set_lock(&vArray[e.source].outEdges.lock);

			//insert outgoing edge
			insertEdge(&vArray[e.source].outEdges, e.destination, e.weight, thInfo[th].edgeCnt);

			//release lock
			omp_unset_lock(&vArray[e.source].outEdges.lock);


			//take a lock on the incoming edge of dest
			omp_set_lock(&vArray[e.destination].inEdges.lock);

			//insert incoming edge
			u64 garbage;
			insertEdge(&vArray[e.destination].inEdges, e.source, e.weight, garbage);

			//release lock
			omp_unset_lock(&vArray[e.destination].inEdges.lock);
		}


//#ifdef _OPENMP
//		int thMask = (1 << getNextPow2Log2(num_threads)) - 1;
//
//		#pragma omp parallel
//		for(u64 i = 0; i < batchSize; i++){
//			const i64 src = el[i].source;
//			const i64 dst = el[i].destination;
//			const i64 actualTh = omp_get_thread_num();
//
//			i64 targetTh = (src / 64) & thMask;
//			if(targetTh == actualTh){
//				if(!el[i].sourceExists){
//					thInfo[actualTh].nodeCnt++;
//				}
//				if(!el[i].destExists){
//					thInfo[actualTh].nodeCnt++;
//				}
//
//				if(!affected[src]){
//					affected[src] = true;
//				}
//
//				//insert out edge
//				insertEdge(&vArray[src].outEdges, dst, el[i].weight, thInfo[actualTh].edgeCnt);
//			}
//
//			targetTh = (dst / 64) & thMask;
//			if(targetTh == actualTh){
//				if(!affected[dst]){
//					affected[dst] = true;
//				}
//				//insert in edge
//				u64 garbase;
//				insertEdge(&vArray[dst].inEdges, src, el[i].weight, garbase);
//			}
//
//		}
//#else
//		for(u64 i = 0; i < batchSize; i++){
//			const u64 src = el[i].source;
//			const u64 dst = el[i].destination;
//
//			if(!el[i].sourceExists){
//				thInfo[0].nodeCnt++;
//			}
//			if(!el[i].destExists){
//				thInfo[0].nodeCnt++;
//			}
//			if(!affected[src]){
//				affected[src] = true;
//			}
//			if(!affected[dst]){
//				affected[dst] = true;
//			}
//
//			//insert out edge
//			insertEdge(&vArray[src].outEdges, dst, el[i].weight, thInfo[0].edgeCnt);
//
//			//insert in edge
//			u64 garbage;
//			insertEdge(&vArray[dst].inEdges, src, el[i].weight, garbage);
//		}
//#endif

		for(u64 i = 0; i < num_threads; i++){
			num_edges += thInfo[i].edgeCnt;
			thInfo[i].edgeCnt = 0;
			num_nodes += thInfo[i].nodeCnt;
			thInfo[i].nodeCnt = 0;
		}
		//num_nodes = el[batchSize - 1].lastAssignedId + 1;
	}
#endif


#ifdef USE_ONLY_HASHMAP

	VertexArray<Neigh> vArray;
	const int num_threads;

	typedef struct{
		u64 edgeCnt = 0;
		u64 nodeCnt = 0;
		u8 pad[48];
	} ThreadInfo;

	alignas(64) ThreadInfo thInfo[32];

	Graphite(bool weighted, bool directed, i64 numNodes, i64 numThreads) : dataStruc(weighted, directed), num_threads(numThreads){
#ifdef _OPENMP
		if(numThreads > 0){
			omp_set_num_threads(numThreads);
		}
#endif

		vArray.resize(numNodes);

		#pragma omp parallel for
		for(u64 i = 0; i < numNodes; i++){
			vArray.outDegree[i] = 0;
			vArray.outCapacity[i] = 0;
			vArray.outNeighArr[i] = nullptr;

			vArray.inDegree[i] = 0;
			vArray.inCapacity[i] = 0;
			vArray.inNeighArr[i] = nullptr;
		}

		//vArray.usingHash.resize(numNodes, false);
		//vArray.dstLocMap.resize(numNodes);

		cout << "Sizeof ThreadInfo: " << sizeof(ThreadInfo) << endl;

		property.resize(numNodes, -1);
		affected.resize(numNodes);
		affected.fill(false);
	}

	~Graphite(){

	}

	void print() override {
		cout << "Done" << endl;
	}

	int64_t in_degree(NodeID n) override {
		return vArray.inDegree[n];
	}

	int64_t out_degree(NodeID n) override {
		return vArray.outDegree[n];
	}


	void insertEdge(u64& deg, u64& cap, Neigh* __restrict &neighs, graphite_hashmap __restrict &locMap, const Idx dstId, const Weight weight, u64& edgeCnt){
		if(deg == cap){
			const u64 newCap = getNextPow2MinRet(cap * 2);

			Neigh* __restrict oldPtr = neighs;
			Neigh* __restrict newPtr = (Neigh*)globalAllocator.allocPow2(newCap * sizeof(Neigh));

			cap = newCap;
			neighs = newPtr;

			if(deg > 0){
				memcpy(newPtr, oldPtr, deg * sizeof(Neigh));
				globalAllocator.freePow2(oldPtr, newCap / 2 * sizeof(Neigh));
			}
		}

		//using hashed mode
		const auto& it = locMap.find(dstId);
		if(it == locMap.end()){
			//edge not found, insert
			Neigh& nn = neighs[deg];
			nn.node = dstId;
			nn.setWeight(weight);
			locMap[dstId] = deg;
			deg++;

			edgeCnt++;
		}
		else{
			//edge found, update weight
			neighs[it->second].setWeight(weight);
		}
	}



	void update(const EdgeList& el) override {
		const u64 batchSize = el.size();

#ifdef _OPENMP
		int thMask = (1 << getNextPow2Log2(num_threads)) - 1;

		#pragma omp parallel
		for(u64 i = 0; i < batchSize; i++){
			const i64 src = el[i].source;
			const i64 dst = el[i].destination;
			const i64 actualTh = omp_get_thread_num();

			i64 targetTh = (src / 64) & thMask;
			if(targetTh == actualTh){
				if(!el[i].sourceExists){
					thInfo[actualTh].nodeCnt++;
				}
				if(!el[i].destExists){
					thInfo[actualTh].nodeCnt++;
				}

				if(!affected[src]){
					affected[src] = true;
				}

				//insert out edge
				insertEdge(vArray.outDegree[src], vArray.outCapacity[src], vArray.outNeighArr[src], vArray.outMapArr[src], dst, el[i].weight, thInfo[actualTh].edgeCnt);
			}

			targetTh = (dst / 64) & thMask;
			if(targetTh == actualTh){
				if(!affected[dst]){
					affected[dst] = true;
				}
				//insert in edge
				u64 garbase;
				insertEdge(vArray.inDegree[dst], vArray.inCapacity[dst], vArray.inNeighArr[dst], vArray.inMapArr[dst], src, el[i].weight, garbase);
			}

		}
#else
		for(u64 i = 0; i < batchSize; i++){
			const u64 src = el[i].source;
			const u64 dst = el[i].destination;

			if(!el[i].sourceExists){
				thInfo[0].nodeCnt++;
			}
			if(!el[i].destExists){
				thInfo[0].nodeCnt++;
			}
			//we do not need atomic operation on affected as long as "some" thread updates it
			if(!affected[src]){
				affected[src] = true;
			}
			if(!affected[dst]){
				affected[dst] = true;
			}

			//insert out edge
			insertEdge(vArray.outDegree[src], vArray.outCapacity[src], vArray.outNeighArr[src], vArray.outMapArr[src], dst, el[i].weight, thInfo[0].edgeCnt);

			//insert in edge
			u64 garbage;
			insertEdge(vArray.inDegree[dst], vArray.inCapacity[dst], vArray.inNeighArr[dst], vArray.inMapArr[dst], src, el[i].weight, garbage);
		}
#endif

		for(u64 i = 0; i < num_threads; i++){
			num_edges += thInfo[i].edgeCnt;
			thInfo[i].edgeCnt = 0;
			num_nodes += thInfo[i].nodeCnt;
			thInfo[i].nodeCnt = 0;
		}
		//num_nodes = el[batchSize - 1].lastAssignedId + 1;
	}


//	void insertEdge(Idx srcId, Idx dstId, Weight weight, int thId){
//		const u64 deg = vArray.outDegree[srcId];
//
//		if(__builtin_expect(deg == vArray.outCapacity[srcId], 0)){
//			const u64 newCap = getNextPow2(vArray.outCapacity[srcId] * 2);
//
//
//			Neigh* __restrict oldPtr = vArray.outNeighArr[srcId];
//			Neigh* __restrict newPtr = (Neigh*)globalAllocator.allocPow2(newCap * sizeof(Neigh));
//
//			vArray.outCapacity[srcId] = newCap;
//			vArray.outNeighArr[srcId] = newPtr;
//
//			if(deg > 0){
//				memcpy(newPtr, oldPtr, deg * sizeof(Neigh));
//				globalAllocator.freePow2(oldPtr, newCap / 2 * sizeof(Neigh));
//			}
//
//			if(__builtin_expect(newCap == HYBRID_HASH_PARTITION, 0)){
//				//switch from linear to hashed mode
//				graphite_hashmap* __restrict hashPtr = (graphite_hashmap*)globalAllocator.allocate(sizeof(graphite_hashmap));
//				new (hashPtr) graphite_hashmap();
//
//				vArray.outMapArr[srcId] = hashPtr;
//
//				//add existing nodes to hash
//				const Neigh* __restrict neighs = newPtr;
//				for(u64 i = 0; i < deg; i++){
//					(*hashPtr)[neighs[i].node] = i;
//				}
//			}
//		}
//
//		//search for existing edge
//		graphite_hashmap* __restrict locMap = vArray.outMapArr[srcId];
//		Neigh* __restrict neigh = vArray.outNeighArr[srcId];
//
//		if(!locMap){
//			//using linear mode
//			Neigh* __restrict nn = nullptr;
//			for(u64 i = 0; i < deg; i++){
//				if(neigh[i].node == dstId){
//					nn = neigh + i;
//					break;
//				}
//			}
//			if(!nn){
//				//edge not found, insert
//				neigh[deg].node = dstId;
//				neigh[deg].setWeight(weight);
//				vArray.outDegree[srcId]++;
//
//				thInfo[thId].edgeCnt++;
//			}
//			else{
//				//edge found, update
//				nn->setWeight(weight);
//			}
//		}
//		else {
//			//using hashed mode
//			const auto& it = locMap->find(dstId);
//			if(it == locMap->end()){
//				//edge not found, insert
//				Neigh& nn = neigh[deg];
//				nn.node = dstId;
//				nn.setWeight(weight);
//				(*locMap)[dstId] = deg;
//				vArray.outDegree[srcId]++;
//
//				thInfo[thId].edgeCnt++;
//			}
//			else{
//				//edge found, update weight
//				neigh[it->second].setWeight(weight);
//			}
//		}
//
//		//we do not need atomic operation on affected as long as "some" thread updates it
//		if(!affected[srcId]){
//			affected[srcId] = true;
//		}
//		if(!affected[dstId]){
//			affected[dstId] = true;
//		}
//	}


#endif


#ifdef USE_ONLY_LINEAR

	VertexArray<Neigh> vArray;
	const int num_threads;

	typedef struct{
		u64 edgeCnt = 0;
		u64 nodeCnt = 0;
		u8 pad[48];
	} ThreadInfo;

	alignas(64) ThreadInfo thInfo[32];

	Graphite(bool weighted, bool directed, i64 numNodes, i64 numThreads) : dataStruc(weighted, directed), num_threads(numThreads){
#ifdef _OPENMP
		if(numThreads > 0){
			omp_set_num_threads(numThreads);
		}
#endif

		vArray.resize(numNodes);

		#pragma omp parallel for
		for(u64 i = 0; i < numNodes; i++){
			vArray.outDegree[i] = 0;
			vArray.outCapacity[i] = 0;
			vArray.outNeighArr[i] = nullptr;

			vArray.inDegree[i] = 0;
			vArray.inCapacity[i] = 0;
			vArray.inNeighArr[i] = nullptr;
		}

		cout << "Sizeof ThreadInfo: " << sizeof(ThreadInfo) << endl;

		property.resize(numNodes, -1);
		affected.resize(numNodes);
		affected.fill(false);
	}

	~Graphite(){

	}

	int64_t in_degree(NodeID n) override {
		return vArray.inDegree[n];
	}

	int64_t out_degree(NodeID n) override {
		return vArray.outDegree[n];
	}


	void insertEdge(u64& deg, u64& cap, Neigh* __restrict &neighs, const Idx dstId, const Weight weight, u64& edgeCnt){
		insTot++;
		if(deg == cap){
			const u64 newCap = getNextPow2MinRet(cap * 2);

			Neigh* __restrict oldPtr = neighs;
			Neigh* __restrict newPtr = (Neigh*)globalAllocator.allocPow2(newCap * sizeof(Neigh));

			cap = newCap;
			neighs = newPtr;

			if(deg > 0){
				memcpy(newPtr, oldPtr, deg * sizeof(Neigh));
				globalAllocator.freePow2(oldPtr, newCap / 2 * sizeof(Neigh));
			}
		}

		//search for existing edge
		Neigh* __restrict nn = nullptr;
		for(u64 i = 0; i < deg; i++){
			if(neighs[i].node == dstId){
				nn = neighs + i;
				break;
			}
		}
		if(!nn){
			//edge not found, insert
			neighs[deg].node = dstId;
			neighs[deg].setWeight(weight);
			deg++;

			edgeCnt++;
			insSucc++;
		}
		else{
			//edge found, update
			nn->setWeight(weight);
		}
	}

	void deleteEdge(u64& deg, u64& cap, Neigh* __restrict &neighs, const Idx dstId, u64& edgeCnt){
		delTot++;
		//search for existing edge
		Neigh* __restrict nn = nullptr;
		for(u64 i = 0; i < deg; i++){
			if(neighs[i].node == dstId){
				nn = neighs + i;
				break;
			}
		}
		if(!nn){
			//edge not found, return
			return;
		}
		else{
			//edge found, delete
			deg--;
			edgeCnt--;
			delSucc++;

			*nn = neighs[deg];
		}

		if((cap > 4) && (deg * 4) <= cap){
			//reduce size
			const u64 newCap = cap / 2;

			Neigh* __restrict oldPtr = neighs;
			Neigh* __restrict newPtr = (Neigh*)globalAllocator.allocPow2(newCap * sizeof(Neigh));

			cap = newCap;
			neighs = newPtr;

			memcpy(newPtr, oldPtr, deg * sizeof(Neigh));
			globalAllocator.freePow2(oldPtr, newCap * 2 * sizeof(Neigh));
		}
	}



	void update(const EdgeList& el) override {
		const u64 batchSize = el.size();

#ifdef _OPENMP
		int thMask = (1 << getNextPow2Log2(num_threads)) - 1;

		#pragma omp parallel
		for(u64 i = 0; i < batchSize; i++){
			const i64 src = el[i].source;
			const i64 dst = el[i].destination;
			const i64 actualTh = omp_get_thread_num();

			i64 targetTh = (src / 64) & thMask;
			if(targetTh == actualTh){
				if(!el[i].sourceExists){
					thInfo[actualTh].nodeCnt++;
				}
				if(!el[i].destExists){
					thInfo[actualTh].nodeCnt++;
				}

				if(!affected[src]){
					affected[src] = true;
				}

				if(!el[i].isDelete){
					//insert out edge
					insertEdge(vArray.outDegree[src], vArray.outCapacity[src], vArray.outNeighArr[src], dst, el[i].weight, thInfo[actualTh].edgeCnt);
				}
				else{
					//delete out edge
					deleteEdge(vArray.outDegree[src], vArray.outCapacity[src], vArray.outNeighArr[src], dst, thInfo[actualTh].edgeCnt);
				}
			}

			targetTh = (dst / 64) & thMask;
			if(targetTh == actualTh){
				if(!affected[dst]){
					affected[dst] = true;
				}

				u64 garbage;

				if(!el[i].isDelete){
					//insert in edge
					insertEdge(vArray.inDegree[dst], vArray.inCapacity[dst], vArray.inNeighArr[dst], src, el[i].weight, garbage);
				}
				else{
					//delete in edge
					deleteEdge(vArray.inDegree[dst], vArray.inCapacity[dst], vArray.inNeighArr[dst], src, garbage);
				}
			}

		}
#else
		for(u64 i = 0; i < batchSize; i++){
			const u64 src = el[i].source;
			const u64 dst = el[i].destination;

			if(!el[i].sourceExists){
				thInfo[0].nodeCnt++;
			}
			if(!el[i].destExists){
				thInfo[0].nodeCnt++;
			}
			//we do not need atomic operation on affected as long as "some" thread updates it
			if(!affected[src]){
				affected[src] = true;
			}
			if(!affected[dst]){
				affected[dst] = true;
			}

			u64 garbage;

			if(!el[i].isDelete){
				//insert out edge
				insertEdge(vArray.outDegree[src], vArray.outCapacity[src], vArray.outNeighArr[src], dst, el[i].weight, thInfo[0].edgeCnt);

				//insert in edge
				insertEdge(vArray.inDegree[dst], vArray.inCapacity[dst], vArray.inNeighArr[dst], src, el[i].weight, garbage);
			}
			else{
				//delete out edge
				deleteEdge(vArray.outDegree[src], vArray.outCapacity[src], vArray.outNeighArr[src], dst, thInfo[0].edgeCnt);

				//delete in edge
				deleteEdge(vArray.inDegree[dst], vArray.inCapacity[dst], vArray.inNeighArr[dst], src, garbage);
			}
		}
#endif

		for(u64 i = 0; i < num_threads; i++){
			num_edges += thInfo[i].edgeCnt;
			thInfo[i].edgeCnt = 0;
			num_nodes += thInfo[i].nodeCnt;
			thInfo[i].nodeCnt = 0;
		}
		//num_nodes = el[batchSize - 1].lastAssignedId + 1;
	}


#endif


#ifdef USE_CAHCE_FRIENDLY_HASH
#include "GraphiteHash.h"

	Vertex<Neigh>* vArray;
	const int num_threads;

	typedef struct{
		u64 edgeCnt = 0;
		u64 nodeCnt = 0;
		u8 pad[48];
	} ThreadInfo;

	alignas(64) ThreadInfo thInfo[32];

	Graphite(bool weighted, bool directed, i64 numNodes, i64 numThreads) : dataStruc(weighted, directed), num_threads(numThreads){
#ifdef _OPENMP
		if(numThreads > 0){
			omp_set_num_threads(numThreads);
		}
#endif
		vArray = (Vertex<Neigh>*)aligned_alloc(64, numNodes * sizeof(Vertex<Neigh>));
		for(u64 i = 0; i < numNodes; i++){
			Vertex<Neigh>& v = vArray[i];
			v.inEdges.degree = 0;
			v.inEdges.capacity = NUM_INITIAL_ELEMS;
			v.inEdges.neighArr = v.inEdges.neigh;

			v.outEdges.degree = 0;
			v.outEdges.capacity = NUM_INITIAL_ELEMS;
			v.outEdges.neighArr = v.outEdges.neigh;
		}

		cout << "Size of GraphiteHash: " << sizeof(GraphiteHash<Neigh>) << endl;
		cout << "Size of Vertex: " << sizeof(Vertex<Neigh>) << endl;
		cout << "Size of ThreadInfo: " << sizeof(ThreadInfo) << endl;
		cout << "Size of Neigh: " << sizeof(Neigh) << endl;

		property.resize(numNodes, -1);
		affected.resize(numNodes);
		affected.fill(false);
	}

	~Graphite(){
		free(vArray);
	}

	void print() override {
		cout << "Done" << endl;
	}

	int64_t in_degree(NodeID n) override {
		return vArray[n].inEdges.degree;
	}

	int64_t out_degree(NodeID n) override {
		return vArray[n].outEdges.degree;
	}

	void update(const EdgeList& el) override {
		const u64 batchSize = el.size();

#ifdef _OPENMP
		int thMask = (1 << getNextPow2Log2(num_threads)) - 1;

		#pragma omp parallel
		for(u64 i = 0; i < batchSize; i++){
			const i64 src = el[i].source;
			const i64 dst = el[i].destination;
			const i64 actualTh = omp_get_thread_num();

			i64 targetTh = (src / 64) & thMask;
			if(targetTh == actualTh){
				if(!el[i].sourceExists){
					thInfo[actualTh].nodeCnt++;
				}
				if(!el[i].destExists){
					thInfo[actualTh].nodeCnt++;
				}

				if(!affected[src]){
					affected[src] = true;
				}

				//insert out edge
				vArray[src].outEdges.insert(dst, thInfo[actualTh].edgeCnt);
			}

			targetTh = (dst / 64) & thMask;
			if(targetTh == actualTh){
				if(!affected[dst]){
					affected[dst] = true;
				}
				//insert in edge
				u64 garbage;
				vArray[dst].inEdges.insert(src, garbage);
			}

		}
#else
		for(u64 i = 0; i < batchSize; i++){
			const u64 src = el[i].source;
			const u64 dst = el[i].destination;

			if(!el[i].sourceExists){
				thInfo[0].nodeCnt++;
			}
			if(!el[i].destExists){
				thInfo[0].nodeCnt++;
			}
			if(!affected[src]){
				affected[src] = true;
			}
			if(!affected[dst]){
				affected[dst] = true;
			}

			//insert out edge
			vArray[src].outEdges.insert(dst, thInfo[0].edgeCnt);

			//insert in edge
			u64 garbage;
			vArray[dst].inEdges.insert(src, garbage);
		}
#endif

		for(u64 i = 0; i < num_threads; i++){
			num_edges += thInfo[i].edgeCnt;
			thInfo[i].edgeCnt = 0;
			num_nodes += thInfo[i].nodeCnt;
			thInfo[i].nodeCnt = 0;
		}
		//num_nodes = el[batchSize - 1].lastAssignedId + 1;
	}

#endif



#ifdef USE_CAHCE_FRIENDLY_HASH_ONLY
#include "GraphiteHash.h"

	Vertex<Neigh>* vArray;
	const int num_threads;

	typedef struct{
		u64 edgeCnt = 0;
		u64 nodeCnt = 0;
		u8 pad[48];
	} ThreadInfo;

	alignas(64) ThreadInfo thInfo[32];

	Graphite(bool weighted, bool directed, i64 numNodes, i64 numThreads) : dataStruc(weighted, directed), num_threads(numThreads){
#ifdef _OPENMP
		if(numThreads > 0){
			omp_set_num_threads(numThreads);
		}
#endif
		vArray = (Vertex<Neigh>*)aligned_alloc(64, numNodes * sizeof(Vertex<Neigh>));
		for(u64 i = 0; i < numNodes; i++){
			Vertex<Neigh>& v = vArray[i];
			v.inEdges.degree = 0;
			v.inEdges.adjList = v.inEdges.neigh;
			//v.inEdges.capacity = NUM_INITIAL_ELEMS * 2;
			//v.inEdges.neighArr = v.inEdges.neigh;

			v.outEdges.degree = 0;
			v.outEdges.adjList = v.outEdges.neigh;
			//v.outEdges.capacity = NUM_INITIAL_ELEMS * 2;
			//v.outEdges.neighArr = v.outEdges.neigh;
		}

		cout << "Size of GraphiteHash: " << sizeof(GraphiteHash<Neigh>) << endl;
		cout << "Size of Vertex: " << sizeof(Vertex<Neigh>) << endl;
		cout << "Size of ThreadInfo: " << sizeof(ThreadInfo) << endl;
		cout << "Size of Neigh: " << sizeof(Neigh) << endl;

		property.resize(numNodes, -1);
		affected.resize(numNodes);
		affected.fill(false);
	}

	~Graphite(){
		free(vArray);
	}

	int64_t in_degree(NodeID n) override {
		return vArray[n].inEdges.degree;
	}

	int64_t out_degree(NodeID n) override {
		return vArray[n].outEdges.degree;
	}

	void update(const EdgeList& el) override {
		const u64 batchSize = el.size();

#ifdef _OPENMP
		int thMask = (1 << getNextPow2Log2(num_threads)) - 1;

		#pragma omp parallel
		for(u64 i = 0; i < batchSize; i++){
			const i64 src = el[i].source;
			const i64 dst = el[i].destination;
			const i64 actualTh = omp_get_thread_num();

			i64 targetTh = (src / 64) & thMask;
			if(targetTh == actualTh){
				if(!el[i].sourceExists){
					thInfo[actualTh].nodeCnt++;
				}
				if(!el[i].destExists){
					thInfo[actualTh].nodeCnt++;
				}

				if(!affected[src]){
					affected[src] = true;
				}

				if(!el[i].isDelete){
					//insert out edge
					vArray[src].outEdges.insert(dst, thInfo[actualTh].edgeCnt);
				}
				else{
					//delete out edge
					vArray[src].outEdges.erase(dst, thInfo[actualTh].edgeCnt);
				}
			}

			targetTh = (dst / 64) & thMask;
			if(targetTh == actualTh){
				if(!affected[dst]){
					affected[dst] = true;
				}

				u64 garbage;
				if(!el[i].isDelete){
					//insert in edge
					vArray[dst].inEdges.insert(src, garbage);
				}
				else{
					//delete in edge
					vArray[dst].inEdges.erase(src, garbage);
				}
			}

		}
#else
		for(u64 i = 0; i < batchSize; i++){
			const u64 src = el[i].source;
			const u64 dst = el[i].destination;

			if(!el[i].sourceExists){
				thInfo[0].nodeCnt++;
			}
			if(!el[i].destExists){
				thInfo[0].nodeCnt++;
			}
			if(!affected[src]){
				affected[src] = true;
			}
			if(!affected[dst]){
				affected[dst] = true;
			}

			u64 garbage;

			if(!el[i].isDelete){
				//insert out edge
				vArray[src].outEdges.insert(dst, thInfo[0].edgeCnt);

				//insert in edge
				vArray[dst].inEdges.insert(src, garbage);
			}
			else{
				//delete out edge
				vArray[src].outEdges.erase(dst, thInfo[0].edgeCnt);

				//delete in edge
				vArray[dst].inEdges.erase(src, garbage);
			}
		}
#endif

		for(u64 i = 0; i < num_threads; i++){
			num_edges += thInfo[i].edgeCnt;
			thInfo[i].edgeCnt = 0;
			num_nodes += thInfo[i].nodeCnt;
			thInfo[i].nodeCnt = 0;
		}
		//num_nodes = el[batchSize - 1].lastAssignedId + 1;
	}

#endif

//	void insertEdge(Idx srcId, Idx dstId, const EProp& eProp){
//		const u64 deg = vArray.getOutDegree(srcId);
//
//		if(__builtin_expect(deg == vArray.getCapacity(srcId), 0)){
//			const u64 newCap = vArray.getCapacity(srcId) * 2;
//
//
//			Idx* 	__restrict dstOld = vArray.getDstIdArray(srcId);
//
//#ifdef HAS_EDGE_PROP
//			EProp* 	__restrict propOld = vSoA.ePropArr[offset];
//			Idx* 	__restrict dstNew = (Idx*)globalAllocator.allocPow2(newCap * (sizeof(Idx) + sizeof(EProp)));
//			EProp* 	__restrict propNew = (EProp*)(dstNew + newCap);
//			//EProp* 	__restrict propNew = (EProp*)globalAllocator.allocPow2(newCap * sizeof(EProp));
//#else
//			Idx* 	__restrict dstNew = (Idx*)globalAllocator.allocPow2(newCap * sizeof(Idx));
//#endif
//
//			for(U64 i = 0; i < deg; i++){
//				dstNew[i] = dstOld[i];
//#ifdef HAS_EDGE_PROP
//				propNew[i] = propOld[i];
//#endif
//			}
//
//#ifdef HAS_EDGE_PROP
//			globalAllocator.freePow2(dstOld, newCap / 2 * (sizeof(Idx) + sizeof(EProp)));
//#else
//			globalAllocator.freePow2(dstOld, newCap / 2 * (sizeof(Idx)));
//#endif
//			//globalAllocator.freePow2(dstOld, newCap / 2 * sizeof(Idx));
//			//globalAllocator.freePow2(propOld, newCap / 2 * sizeof(EProp));
//
//#if defined(USE_FULL_HASHMAP) || defined(USE_HYBRID_HASH_V2)
//			vArray.dstLocMap[srcId][0].reserve(newCap);
//#endif
//			vArray.setCapacity(srcId, newCap);
//			vArray.setDstIdArray(srcId, dstNew);
//#ifdef HAS_EDGE_PROP
//			vSoA.ePropArr[offset] = propNew;
//#endif
//		}
//
//#if defined(USE_FULL_HASHMAP)
//		vSoA.dstLocMap[offset][dstId] = deg;
//#elif defined(USE_HYBRID_HASHMAP)
//		if(deg >= HYBRID_HASH_PARTITION) {
//			vSoA.dstLocMap[offset][dstId] = deg;
//		}
//#elif defined(USE_HYBRID_HASH_V2)
//		if(vArray.dstLocMap[srcId] != nullptr) {
//			vArray.dstLocMap[srcId][0][dstId] = deg;
//		}
//		else if(deg >= HYBRID_HASH_THRESHOLD){
//			vArray.dstLocMap[srcId] = new graphite_hashmap;
//			vArray.dstLocMap[srcId][0][dstId] = deg;
//		}
//
//#endif
//
//		vArray.getDstIdArray(srcId)[deg] = dstId;
//#ifdef HAS_EDGE_PROP
//		vSoA.ePropArr[offset][deg] = eProp;
//#endif
//		vArray.getOutDegree(srcId)++;
//	}



//	void deleteEdge(Idx srcId, Idx dstId){
//		Idx* __restrict dstIdArr = vArray.getDstIdArray(srcId);
//
//#ifdef HAS_EDGE_PROP
//		EProp* __restrict ePropArr = vSoA.ePropArr[offset];
//#endif
//
//		if(vArray.getOutDegree(srcId)){
//#if defined(USE_FULL_HASHMAP)
//			const auto& dstIt = vSoA.dstLocMap[offset].find(dstId);
//			if(dstIt != vSoA.dstLocMap[offset].end()){
//				const Idx dstLoc = dstIt->second;
//				vSoA.dstLocMap[offset].erase(dstIt);
//				const Idx outDeg = --vSoA.outDegree[offset];
//				if(__builtin_expect(outDeg != dstLoc, 1)){		//not the last element
//					dstIdArr[dstLoc] = dstIdArr[outDeg];
//#ifdef HAS_EDGE_PROP
//					ePropArr[dstLoc] = ePropArr[outDeg];
//#endif
//					vSoA.dstLocMap[offset][dstIdArr[outDeg]] = dstLoc;
//				}
//			}
//#elif defined(USE_SORTED_EDGES)
//			// First thing first, find the location of the dst node
//			const Idx outDeg = vSoA.outDegree[offset] - 1;
//			const Idx lastDstId = dstIdArr[outDeg];
//			if(__builtin_expect(lastDstId != dstId, 1)){		//not the last element
//				const u64 sortedLen = vSoA.sortedLen[offset];
//				u64 dstLoc = lower_bound(dstIdArr, dstIdArr + sortedLen, dstId) - dstIdArr;
//				if((dstLoc == sortedLen) || (dstIdArr[dstLoc] != dstId)){		//not found in the sorted region, linear search rest
//					dstLoc = sortedLen - 1;
//					while(dstIdArr[++dstLoc] != dstId && dstLoc < outDeg);
//				}
//
//				if(__builtin_expect(dstLoc != outDeg, 1)){
//					//found
//					dstIdArr[dstLoc] = dstIdArr[outDeg];
//#ifdef HAS_EDGE_PROP
//					ePropArr[dstLoc] = ePropArr[outDeg];
//#endif
//					vSoA.outDegree[offset]--;
//				}
//			}
//			else{ //last element
//				vSoA.outDegree[offset]--;
//			}
//
//			// sorted len goes over the out degree. Adjust.
//			if(vSoA.sortedLen[offset] > vSoA.outDegree[offset]){
//				vSoA.sortedLen[offset]--;
//			}
//
//#elif defined(USE_HYBRID_HASHMAP)
//			const u64 outDeg = vSoA.outDegree[offset] - 1;
//
//			if(dstIdArr[outDeg] == dstId){	//last element
//				vSoA.outDegree[offset]--;
//				if(outDeg >= HYBRID_HASH_PARTITION){
//					vSoA.dstLocMap[offset].erase(dstId);
//				}
//				return;
//			}
//
//
//			u64 dstLoc = -1;
//
//			if(outDeg >= HYBRID_HASH_PARTITION){
//				const auto& dstIt = vSoA.dstLocMap[offset].find(dstId);
//				if(dstIt != vSoA.dstLocMap[offset].end()){
//					//found
//					dstLoc = dstIt->second;
//					vSoA.dstLocMap[offset].erase(dstIt);
//				}
//			}
//
//			if(dstLoc == -1){ //not yet found, do linear search
//				while(dstIdArr[++dstLoc] != dstId && dstLoc < outDeg);
//				if(__builtin_expect(dstLoc == outDeg, 0)){
//					//not found
//					return;
//				}
//			}
//
//
//
//			//valid dstLoc
//			dstIdArr[dstLoc] = dstIdArr[outDeg];
//#ifdef HAS_EDGE_PROP
//			ePropArr[dstLoc] = ePropArr[outDeg];
//#endif
//			vSoA.dstLocMap[offset][dstIdArr[outDeg]] = dstLoc;
//
//#elif defined(USE_HYBRID_HASH_V2)
//			if(vArray.dstLocMap[srcId] == nullptr){
//				//Not using hashmap
//				const Idx outDeg = vArray.getOutDegree(srcId) - 1;
//				Idx lastDstId = dstIdArr[outDeg];
//				if(__builtin_expect(lastDstId != dstId, 1)){		//not the last element
//					u64 dstLoc = -1;
//					while(dstIdArr[++dstLoc] != dstId && dstLoc < outDeg);
//					if(__builtin_expect(dstLoc != outDeg, 1)){
//						dstIdArr[dstLoc] = dstIdArr[outDeg];
//#ifdef HAS_EDGE_PROP
//						ePropArr[dstLoc] = ePropArr[outDeg];
//#endif
//						vArray.getOutDegree(srcId)--;
//					}
//				}
//				else{
//					vArray.getOutDegree(srcId)--;		//last element
//				}
//			}
//			else{
//				//Using hashmap
//				const auto& dstIt = vArray.dstLocMap[srcId]->find(dstId);
//				if(dstIt != vArray.dstLocMap[srcId]->end()){
//					const Idx dstLoc = dstIt->second;
//					vArray.dstLocMap[srcId]->erase(dstIt);
//					const Idx outDeg = --vArray.outDegree[srcId];
//					if(__builtin_expect(outDeg != dstLoc, 1)){		//not the last element
//						dstIdArr[dstLoc] = dstIdArr[outDeg];
//#ifdef HAS_EDGE_PROP
//						ePropArr[dstLoc] = ePropArr[outDeg];
//#endif
//						vArray.dstLocMap[srcId][0][dstIdArr[outDeg]] = dstLoc;
//					}
//				}
//			}
//
//#else
//			const Idx outDeg = vArray.getOutDegree(srcId) - 1;
//			Idx lastDstId = dstIdArr[outDeg];
//			if(__builtin_expect(lastDstId != dstId, 1)){		//not the last element
//				u64 dstLoc = -1;
//				while(dstIdArr[++dstLoc] != dstId && dstLoc < outDeg);
//				if(__builtin_expect(dstLoc != outDeg, 1)){
//					dstIdArr[dstLoc] = dstIdArr[outDeg];
//#ifdef HAS_EDGE_PROP
//					ePropArr[dstLoc] = ePropArr[outDeg];
//#endif
//					vArray.getOutDegree(srcId)--;
//				}
//			}
//			else{
//				vArray.getOutDegree(srcId)--;		//last element
//			}
//#endif
//		}
//	}

//
//	void insertEdges(Idx srcId, U64 count, const Idx* __restrict dstIds, const EProp* __restrict eProps){
//		Vertex& vv = vertices[srcId];
//		VExtra& ve = vExtras[srcId];
//
//		if(ve.eCap < (vv.outDeg + count)){
//			ve.eCap = getNextPow2(vv.outDeg + count);
//			Idx* 	__restrict dstOld = vv.dstArr;
//			EProp* 	__restrict propOld = vv.ePropArr;
//			vv.dstArr = (Idx*)aligned_alloc(64, ve.eCap * sizeof(Idx));
//			vv.ePropArr = (EProp*)aligned_alloc(64, ve.eCap * sizeof(EProp));
//			#pragma omp parallel for simd
//			for(U64 i = 0; i < vv.outDeg; i++){
//				vv.dstArr[i] = dstOld[i];
//				vv.ePropArr[i] = propOld[i];
//			}
//			free(dstOld);
//			free(propOld);
//			ve.dstLocMap.reserve(ve.eCap);
//		}
//		#pragma omp parallel for
//		for(U64 i = 0; i < count; i++){
//			vv.ePropArr[vv.outDeg + i] = eProps[i];
//			vv.dstArr[vv.outDeg + i] = dstIds[i];
//			ve.dstLocMap[dstIds[i]] = i;
//		}
//		vv.outDeg += count;
//	}


	//typedef const vProp_t& (*processEdgeFunc)(const eProp_t& eProp, const vProp_t& srcProp, const vProp_t& dstProp);
	//typedef const vProp_t& (*reduceFunc)(const vProp_t& temp, const vProp_t& res);

	//---------------------------------- VERTEX ------------------------------------------
	/*inline I64 insertVertex(const VProp& vProp = 0){
		if(deletedVertices.empty()){
			vProps.push_back(vProp);
			outEdges.push_back(new EdgeArrayBlock<eProp_t>);
			//valid.push_back(true);
			return (vProps.size() - 1);
		}
		else{
			I64 vidx = deletedVertices.top();
			//valid[vidx] = true;
			deletedVertices.pop();
			return vidx;
		}
	}*/


	/*U64 insertManyVertex(U64 count){
		U64 startId = vertices.size();
		vertices.resize(startId + count);
		for(U64 i = 0; i < count; i++){
			outEdges.push_back(new EdgeArrayBlock<eProp_t>);
			//outEdges.push_back(new EdgeArrayVector<eProp_t>);
		}
		return startId;
	}*/


//	U64 insertManyVertex(const U64 count, const VProp& vProp = 0){
//		const U64 startId = vSize;
//		const U64 endId = vSize + count;
//		if(vCap < endId){
//			vCap = getNextPow2(endId);
//			Vertex* __restrict newVer = (Vertex*)aligned_alloc(64, vCap * sizeof(Vertex));
//			VExtra* __restrict newExt = (VExtra*)aligned_alloc(64, vCap * sizeof(VExtra));
//			#pragma omp parallel for simd
//			for(U64 i = 0; i < vSize; i++){
//				newVer[i] = vertices[i];
//				newExt[i] = vExtras[i];
//			}
//			free(vertices);
//			free(vExtras);
//			vertices = newVer;
//			vExtras = newExt;
//		}
//		#pragma omp parallel for
//		for(U64 i = startId; i < endId; i++){
//			vertices[i].prop = vProp;
//			vertices[i].outDeg = 0;
//			vExtras[i].eCap = getNextPow2(0);
//			vertices[i].dstArr = (Idx*)aligned_alloc(64, vExtras[i].eCap * sizeof(Idx));
//			vertices[i].ePropArr = (EProp*)aligned_alloc(64, vExtras[i].eCap * sizeof(EProp));
//			new (&vExtras[i].dstLocMap) std::unordered_map<Idx, U64>;
//		}
//		vSize = endId;
//		return startId;
//	}

	/*inline void deleteVertex(I64 vIdx){
		//TODO
		const I64 cnt = valid.size();
		for(I64 i = 0; i < cnt; i++){
			if(valid[i]){
				outEdges[i]->deleteEdge(vIdx);
			}
		}
	}*/

	/*inline vProp_t& getVProp(I64 vIdx) {
		return vProps[vIdx];
	}

	inline const vProp_t& getVProp(I64 vIdx) const {
		return vProps[vIdx];
	}

	inline void setVProp(I64 vIdx, const vProp_t &vProp){
		vProps[vIdx] = vProp;
	}

	inline U64 getVCount() const {
		return (vProps.size() - deletedVertices.size());
	}

	//----------------------------------- EDGE -------------------------------------------
	inline void insertEdge(I64 srcId, I64 dstId, const eProp_t &eProp = {}){
		//eCount++;
		outEdges[srcId]->insertEdge(dstId, eProp);
	}

	inline void updateEdge(I64 srcId, I64 dstId, const eProp_t &eProp){
		outEdges[srcId]->updateEdge(dstId, eProp);
	}

	inline void deleteEdge(I64 srcId, I64 dstId){
		outEdges[srcId]->deleteEdge(dstId);
		//eCount--;
	}

	inline eProp_t& getEProp(I64 srcId, I64 dstId){
		return outEdges[srcId]->getEProp(dstId);
	}

	inline const eProp_t& getEProp(I64 srcId, I64 dstId) const {
		return outEdges[srcId]->getEProp(dstId);
	}

	inline U64 getECount() const {
		return eCount;
	}

	inline auto getEdgeArray(I64 srcId){
		return outEdges[srcId];
	}*/

	//----------------------------------- ANALYTICS --------------------------------------
	/*void run_analysis_algo(processEdgeFunc processEdge, reduceFunc reduce){
		while(!activeVertices.empty()){
			processing_phase(processEdge, reduce);
			commit_phase();
		}
	}*/

//	void run_bfs(I64 rootIdx){
//		typedef struct {
//			U64 	outDeg;
//			Idx* 	dstArr;
//		} EdgeInfo;
//
//		const U64 maxTh = omp_get_max_threads();
//
//		std::vector<EdgeInfo> q;
//		U64 qTail = 0;
//
//		q.push_back({vertices[rootIdx].outDeg, vertices[rootIdx].dstArr});
//		vertices[rootIdx].prop = 0;
//
//		U64 iter = 1;
//		while((q.size() - qTail) < maxTh) {
//			const U64 qHeadOrig = q.size();
//			for(U64 qT = qTail; qT < qHeadOrig; qT++){
//				//pop qTail
//				const EdgeInfo oe = q[qT];
//				for(U64 i = 0; i < oe.outDeg; i++){
//					const I64 dstId = oe.dstArr[i];
//					if(vertices[dstId].prop > iter){
//						vertices[dstId].prop = iter;
//						q.push_back({vertices[dstId].outDeg, vertices[dstId].dstArr});
//					}
//				}
//			}
//			qTail = qHeadOrig;
//			iter++;
//		}
//
//		typedef struct {
//			std::vector<EdgeInfo>	q;
//			U64						qTail;
//			U64						iter;
//			U8						pad[24];
//		} ThQ;
//
//		std::vector<ThQ> qVec(maxTh);
//
//		// Initialize
//		for(U64 i = 0; i < maxTh; i++){
//			qVec[i].qTail = 0;
//			qVec[i].iter = iter;
//		}
//
//		//distribute vertices
//		U64 idx = 0;
//		for(U64 i = qTail; i < q.size(); i++){
//			if(idx == maxTh){
//				idx = 0;
//			}
//			qVec[idx].q.push_back(q[i]);
//			idx++;
//		}
//
//
//		#pragma omp parallel
//		{
//			U64 th = omp_get_thread_num();
//			std::vector<EdgeInfo>& q = qVec[th].q;
//			U64&	qTail = qVec[th].qTail;
//			U64&	iter = qVec[th].iter;
//			while(q.size() - qTail) {
//				const U64 qHeadOrig = q.size();
//				for(U64 qT = qTail; qT < qHeadOrig; qT++){
//					//pop qTail
//					const EdgeInfo oe = q[qT];
//					for(U64 i = 0; i < oe.outDeg; i++){
//						const I64 dstId = oe.dstArr[i];
//						if(vertices[dstId].prop > iter){
//							vertices[dstId].prop = iter;
//							q.push_back({vertices[dstId].outDeg, vertices[dstId].dstArr});
//						}
//					}
//				}
//				qTail = qHeadOrig;
//				iter++;
//			}
//		}
//	}
//
//	void run_bfs_old(I64 rootIdx){
//		typedef struct {
//			U64 	outDeg;
//			Idx* 	dstArr;
//		} EdgeInfo;
//
//		EdgeInfo* __restrict q = (EdgeInfo*)aligned_alloc(64, vSize * sizeof(EdgeInfo));
//		U64 qHead = 1;
//		U64 qTail = 0;
//
//		q[0] = {vertices[rootIdx].outDeg, vertices[rootIdx].dstArr};
//		vertices[rootIdx].prop = 0;
//
//		U64 iter = 1;
//		while(qHead - qTail) {
//			const U64 qHeadOrig = qHead;
//			#pragma omp parallel for
//			for(U64 qT = qTail; qT < qHeadOrig; qT++){
//				//pop qTail
//				const EdgeInfo& oe = q[qT];
//				for(U64 i = 0; i < oe.outDeg; i++){
//					const I64 dstId = oe.dstArr[i];
//					if(vertices[dstId].prop > iter){
//						vertices[dstId].prop = iter;
//						U64 loc;
//						#pragma omp atomic capture
//						loc = qHead++;
//						q[loc] = {vertices[dstId].outDeg, vertices[dstId].dstArr};
//					}
//				}
//			}
//			qTail = qHeadOrig;
//			iter++;
//		}
//	}

	/*void addVertexToWavefront(I64 vIdx){
		activeVertices.push_back(vIdx);
	}

	inline void processing_phase(processEdgeFunc processEdge, reduceFunc reduce){
		const U64 activeCount = activeVertices.size();
		#pragma omp parallel for
		for(U64 i = 0; i < activeCount; i++){
			I64 srcId = activeVertices[i];
			auto outEdgeArr = outEdges[srcId];
			const U64 outEdgeCount = outEdgeArr.size();
			const vProp_t& srcProp = vProps[srcId];
			for(U64 i = 0; i < outEdgeCount; i++){
				const auto edge = outEdgeArr->edges[i];
				const vProp_t& res = processEdge(edge.eProp, srcProp, vProps[edge.dstId]);
				vTempProps[edge.dstId] = reduce(vTempProps[edge.dstId], res);
			}
		}
	}

	inline void commit_phase(){
		const U64 vertexCount = getVCount();
		#pragma omp parallel for
		for(U64 i = 0; i < vertexCount; i++){
			if(vProps[i] != vTempProps[i]){
				vProps[i] = vTempProps[i];
				omp_set_lock(&lock);
				activeVertices.push_back(i);
				omp_unset_lock(&lock);
			}
		}
	}*/

	//----------------------------------- OTHERS -----------------------------------------
	/*void compact(){

	}

	void garbage_collection(){

	}*/

	//static Graphite* buildFromTextCSR(const std::string &fname){

	//}



	/*static u64 binSearch(const i64* sortedArr, i64 val, u64 s, u64 e){
		if((e - s) <= 1){
			i64 diffE = sortedArr[e * VECTOR_WIDTH] - val;
			i64 diffS = val - sortedArr[s * VECTOR_WIDTH];
			if(diffE < diffS){
				return e;
			}
			return s;
		}
		u64 mid = (s + e) / 2;
		u64 midVal = sortedArr[mid * VECTOR_WIDTH];
		if(val == midVal){
			return mid;
		}
		else if(val < midVal){
			return binSearch(sortedArr, val, s, mid);
		}
		else {
			return binSearch(sortedArr, val, mid, e);
		}
	}*/

//	static Graphite* buildFromStingerCSR(
//			const i64 numThreads,
//			const U64 nv,
//			const U64 ne,
//			const I64* __restrict off,
//			const I64* __restrict ind,
//			const I64* __restrict weight,
//			const VProp& defaultProp = 0)
//	{
//		auto gra = new Graphite<VProp, EProp>(numThreads);
//
//		gra->vArray.resize(nv);
//
//
////
////		ofstream outf("edgelist.txt");
////		for(u64 v = 0; v < nv; v++){
////			const U64 deg = off[v + 1] - off[v];
////			const I64* __restrict indE = ind + off[v];
////			const I64* __restrict wgtE = weight + off[v];
////			for(u64 e = 0; e < deg; e++){
////				outf << v << " " << indE[e] << "\n";
////			}
////		}
////		outf.close();
//
////		struct {
////			u64 e = 0;
////			u8 pad[56];
////		} cnt[25];
////
////
////		const u64 perThEdge = ne / omp_get_max_threads();
////		const u64 nvVect = (nv + VECTOR_WIDTH - 1)/VECTOR_WIDTH;
////		#pragma omp parallel
////		{
////			u64 th = omp_get_thread_num();
////			u64 target = (th + 1) * perThEdge;
////			cnt[th+1].e = binSearch(off, target, 0, nvVect) * VECTOR_WIDTH;
////		}
////
////		auto& vArr = gra->vertexArray;
////		#pragma omp parallel
////		{
////			const u64 th = omp_get_thread_num();
////			for(u64 v = cnt[th].e; v < cnt[th+1].e; v++){
////				VertexSoA<VProp, EProp>& vSoA = vArr[v / VECTOR_WIDTH];
////				const u64 offset = v % VECTOR_WIDTH;
////				const U64 deg = off[v + 1] - off[v];
////				vSoA.outDegree[offset] = deg;
////
////				const u64 cap = std::max(getNextPow2(deg), 4UL);
////				vSoA.capacity[offset] = cap;
////
////				vSoA.dstIdArr[offset] = (Idx*)globalAllocator.allocPow2(cap * sizeof(Idx));
////				vSoA.ePropArr[offset] = (EProp*)globalAllocator.allocPow2(cap * sizeof(EProp));
////				//vArr[v / VECTOR_WIDTH].dstIdArr[v % VECTOR_WIDTH] = (Idx*)aligned_alloc(64, cap * sizeof(Idx));
////				//vArr[v / VECTOR_WIDTH].ePropArr[v % VECTOR_WIDTH] = (EProp*)aligned_alloc(64, cap * sizeof(EProp));
////
////				new (&vSoA.dstLocMap[offset]) unordered_map<Idx, U64>;
////				vSoA.dstLocMap[offset].reserve(deg);
////
////				const I64* __restrict indE = ind + off[v];
////				const I64* __restrict wgtE = weight + off[v];
////
////				Idx* dstIdArr = vSoA.dstIdArr[offset];
////				EProp* ePropArr = vSoA.ePropArr[offset];
////				auto& dstLocMap = vSoA.dstLocMap[offset];
////
////				for(u64 e = 0; e < deg; e++){
////					ePropArr[e] = wgtE[e];
////					dstIdArr[e] = indE[e];
////					dstLocMap[ind[e]] = e;
////				}
////
////			}
////		}
//
//		auto& vArr = gra->vArray;
//		#pragma omp parallel for
//		for(u64 v = 0; v < nv; v++){
//			const U64 deg = off[v + 1] - off[v];
//			vArr.setOutDegree(v, deg);
//			vArr.setVProp(v, defaultProp);
//			const u64 cap = getNextPow2(deg);
//			vArr.setCapacity(v, cap);
//
//#ifdef HAS_EDGE_PROP
//			vSoA.dstIdArr[offset] = (Idx*)globalAllocator.allocPow2(cap * (sizeof(Idx) + sizeof(EProp)));
//			vSoA.ePropArr[offset] = (EProp*)(vSoA.dstIdArr[offset] + cap);
//			//vSoA.ePropArr[offset] = (EProp*)globalAllocator.allocPow2(cap * sizeof(EProp));
//#else
//			//vSoA.dstIdArr[offset] = (Idx*)globalAllocator.allocPow2(cap * sizeof(Idx));
//			vArr.setDstIdArray(v, (Idx*)globalAllocator.allocPow2(cap * sizeof(Idx)));
//#endif
//
//			const I64* __restrict indE = ind + off[v];
//			const I64* __restrict wgtE = weight + off[v];
//
//			I64* __restrict dstIdArr = vArr.getDstIdArray(v);
//#ifdef HAS_EDGE_PROP
//			I64* __restrict ePropArr = vSoA.ePropArr[offset];
//#endif
//
//
//#if defined(USE_FULL_HASHMAP)
//			auto& dstLocMap = vSoA.dstLocMap[offset];
//			new (&dstLocMap) graphite_hashmap;
//			dstLocMap.reserve(cap);
//			for(u64 e = 0; e < deg; e++){
//#ifdef HAS_EDGE_PROP
//				ePropArr[e] = wgtE[e];
//#endif
//				dstIdArr[e] = indE[e];
//				dstLocMap[ind[e]] = e;
//			}
//#elif defined(USE_HYBRID_HASHMAP)
//			u64 linearRegion = std::min(deg, HYBRID_HASH_PARTITION);
//			for(u64 e = 0; e < linearRegion; e++){
//#ifdef HAS_EDGE_PROP
//				ePropArr[e] = wgtE[e];
//#endif
//				dstIdArr[e] = indE[e];
//			}
//			auto& dstLocMap = vSoA.dstLocMap[offset];
//			new (&dstLocMap) graphite_hashmap;
//			if(deg > HYBRID_HASH_PARTITION){
//				dstLocMap.reserve(deg - HYBRID_HASH_PARTITION);
//				for(u64 e = HYBRID_HASH_PARTITION; e < deg; e++){
//#ifdef HAS_EDGE_PROP
//					ePropArr[e] = wgtE[e];
//#endif
//					dstIdArr[e] = indE[e];
//					dstLocMap[ind[e]] = e;
//				}
//			}
//#elif defined(USE_SORTED_EDGES)
//			for(u64 e = 0; e < deg; e++){
//#ifdef HAS_EDGE_PROP
//				ePropArr[e] = wgtE[e];
//#endif
//				dstIdArr[e] = indE[e];
//			}
//			sort(dstIdArr, dstIdArr + deg);
//			vSoA.sortedLen[offset] = deg;
//#elif defined(USE_HYBRID_HASH_V2)
//			if(deg < HYBRID_HASH_THRESHOLD){
//				//Do not use hashmap
//				vArr.dstLocMap[v] = nullptr;
//				for(u64 e = 0; e < deg; e++){
//					const Idx dstId = indE[e];
//#ifdef HAS_EDGE_PROP
//					ePropArr[e] = wgtE[e];
//#endif
//					dstIdArr[e] = dstId;
//				}
//			}
//			else{
//				//Use hashmap
//				//vArr.dstLocMap[v] = (graphite_hashmap*)globalAllocator.allocate(sizeof(graphite_hashmap));
//				//new (&vArr.dstLocMap[v]) graphite_hashmap();
//				vArr.dstLocMap[v] = new graphite_hashmap();
//				for(u64 e = 0; e < deg; e++){
//					const Idx dstId = indE[e];
//#ifdef HAS_EDGE_PROP
//					ePropArr[e] = wgtE[e];
//#endif
//					dstIdArr[e] = dstId;
//					vArr.dstLocMap[v][0][dstId] = e;
//					//vArr.dstLocMap[v]->at(dstId) = e;
//				}
//			}
//
//#else
//			for(u64 e = 0; e < deg; e++){
//				const Idx dstId = indE[e];
//#ifdef HAS_EDGE_PROP
//				ePropArr[e] = wgtE[e];
//#endif
//				dstIdArr[e] = dstId;
//#ifdef USE_INCOMING_EDGES
//				omp_set_lock(&vArr.lock[dstId]);
//				vArr.inEdges[dstId].push_back(v);
//				omp_unset_lock(&vArr.lock[dstId]);
//#endif
//			}
//#if defined(SORT_EDGES_AT_BUILD)
//			std::sort(dstIdArr, dstIdArr + deg);
//#endif
//#endif
//		}
//
//		return gra;
//	}

	/*static Graphite* buildFromCSC(const std::string &fname){

	}

	static Graphite* buildFromCOO(const std::string &fname){

	}*/


//	void streamStingerActions(U64 naction, const I64* __restrict__ actions){
//
//		std::vector<std::pair<u64, u64> > srcDst;
//		for(u64 i = 0; i < naction; i++){
//			long src = actions[i*2];
//			if(src >= 0){
//				srcDst.push_back({actions[i*2], actions[i*2 + 1]});
//			}
//			else{
//				srcDst.push_back({~actions[i*2], ~actions[i*2 + 1]});
//			}
//
//		}
//
//		cout << srcDst.size() << endl;
//		cout << "===========insertions============" << endl;
//
//		u64 iter = 0;
//		u64 totEdge = 0;
//		double totElp = 0.0;
//
//		u64 nv = vArray.vSize;
//
//		while(1){
//			iter++;
//			if(iter == 50){
//				break;
//			}
//
//			tic();
//			for(u64 i = 0; i < 1024*1024; i++){
//				insertEdge(srcDst[totEdge].first, srcDst[totEdge].second, 1);
//				totEdge++;
//			}
//			totElp += toc();
//
//			cout << totEdge / 1000000.0 / totElp << endl;
//		}
//
//
//
//
//		cout << "===========Deletions============" << endl;
//
//		iter = 0;
//		totEdge = 0;
//		totElp = 0.0;
//
//		while(1){
//			iter++;
//			if(iter == 50){
//				break;
//			}
//
//			tic();
//			for(u64 i = 0; i < 1024*1024; i++){
//				deleteEdge(srcDst[totEdge].first, srcDst[totEdge].second);
//				totEdge++;
//			}
//			totElp += toc();
//
//			cout << totEdge / 1000000.0 / totElp << endl;
//		}


//		//Divide work
//		int numThreads = omp_get_max_threads();
//		int maxTh = 1;
//		while(maxTh <= numThreads){
//			maxTh *= 2;
//		}
//		maxTh--;
//
//		#pragma omp parallel num_threads(maxTh + 1)
//		for(U64 i = 0; i < naction; i++) {
//			I64 src = actions[i * 2];
//			const I64 mask = 0UL - ((src >> 63) & 1UL);
//			src = src ^ mask;
//			int mappedTh = src & maxTh;
//			if(omp_get_thread_num() == mappedTh){
//				I64 dst = actions[i * 2 + 1] ^ mask;
//				if(mask == 0UL){
//					insertEdge(src, dst, 1);
//				}
//				else{
//					deleteEdge(src, dst);
//				}
//			}
//		}
//	}

	//similar to what is done by stinger (queue based)
//	u64 run_bfs_v1(const Idx& root){
//		Idx* __restrict activeList = (Idx*)globalAllocator.allocate(vArray.vSize * sizeof(Idx));
//		u64 qHead = 1;
//		u64 qTail = 0;
//
//		activeList[0] = root;
//		vArray.setVProp(root, 0);
//
//		U64 iter = 1;
//		while(qHead - qTail) {
//			const U64 qHeadOrig = qHead;
//			#pragma omp parallel for
//			for(U64 qT = qTail; qT < qHeadOrig; qT++){
//				//pop qTail
//				const u64 vIdx = activeList[qT];
//				const u64 outDeg = vArray.outDegree[vIdx];
//				const Idx* __restrict dstIdArr = vArray.dstIdArr[vIdx];
//				for(U64 i = 0; i < outDeg; i++){
//					const I64 dstId = dstIdArr[i];
//					if(vArray.prop[dstId] > iter){
//						vArray.prop[dstId] = iter;
//						U64 loc;
//						#pragma omp atomic capture
//						loc = qHead++;
//						activeList[loc] = dstId;
//					}
//				}
//			}
//			qTail = qHeadOrig;
//			iter++;
//		}
//		globalAllocator.deallocate(activeList, vArray.vSize * sizeof(Idx));
//		return iter;
//	}
//
//
//	u64 run_bfs_v2(const Idx& root){
//		const u64 numV = vArray.vSize;
//
//		vector<bool> currAct(numV, 0);
//		vector<bool> nextAct(numV, 0);
//
//
//		vArray.setVProp(root, 0);
//		currAct[root] = 1;
//
//		u64 iter = 1;
//		bool isChanged = true;
//		while(isChanged){
//			isChanged = false;
//			#pragma omp parallel for num_threads(64)
//			for(u64 v = 0; v < vArray.vSize; v++){
//				if(currAct[v]){		//active
//					currAct[v] = 0;
//					const u64 outDeg = vArray.outDegree[v];
//					const Idx* __restrict dstIdArr = vArray.dstIdArr[v];
//					//how about using AVX here?
//					for(U64 i = 0; i < outDeg; i++){
//						const I64 dstId = dstIdArr[i];
//						if(vArray.prop[dstId] > iter){
//							vArray.prop[dstId] = iter;
//							// make dstId active for next iter
//							if(nextAct[dstId] == 0){	//checking this condition reduces thrashing
//								nextAct[dstId] = 1;
//							}
//							isChanged = true;
//						}
//					}
//				}
//			}
//			swap(currAct, nextAct);
//			iter++;
//		}
//
//		return iter;
//	}
//
//
//	u64 run_bfs_v3(const Idx& root){
//		const u64 numV = vArray.vSize;
//
//		vector<bool> currAct(numV, 0);		// wavefront for the current iter
//		vector<bool> nextAct(numV, 0);		// wavefront for the next iter
//
//		currAct[root] = 1;
//
//		u64 iter = 0;
//		bool isChanged = true;
//		while(isChanged){
//			isChanged = false;
//			#pragma omp parallel for num_threads(32)
//			for(u64 v = 0; v < vArray.vSize; v++){
//				if(__builtin_expect(currAct[v], 0)){
//					currAct[v] = 0;
//					if(vArray.prop[v] > iter){
//						vArray.prop[v] = iter;
//						isChanged = true;
//						const u64 outDeg = vArray.outDegree[v];
//						const Idx* __restrict dstIdArr = vArray.dstIdArr[v];
//						for(U64 i = 0; i < outDeg; i++){
//							if(nextAct[dstIdArr[i]] == 0){		//checking this condition should reduces thrashing
//								nextAct[dstIdArr[i]] = 1;
//							}
//						}
//					}
//				}
//			}
//			swap(currAct, nextAct);
//			iter++;
//		}
//
//		return iter;
//	}
//
//
//	u64 run_bfs_v4(const Idx& root){
//		const u64 numV = vArray.vSize;
//
//		vector<bool> currAct(numV, 0);		// wavefront for the current iter
//		vector<bool> nextAct(numV, 0);		// wavefront for the next iter
//		vector<bool> visited(numV, 0);		// records visited vertices
//
//		currAct[root] = 1;
//
//		u64 iter = 0;
//		bool isChanged = true;
//		while(isChanged){
//			isChanged = false;
//			#pragma omp parallel for num_threads(32)
//			for(u64 v = 0; v < vArray.vSize; v++){
//				if (currAct[v]){
//					currAct[v] = 0;
//					vArray.prop[v] = iter;	//update value
//					isChanged = true;
//					const u64 outDeg = vArray.outDegree[v];
//					const Idx* __restrict dstIdArr = vArray.dstIdArr[v];
//					for(U64 e = 0; e < outDeg; e++){
//						const u64 dstId = dstIdArr[e];
//						if(visited[dstId] == 0){
//							visited[dstId] = 1;
//							nextAct[dstId] = 1;
//						}
//					}
//				}
//			}
//			swap(currAct, nextAct);
//			iter++;
//		}
//
//		return iter;
//	}
//
//
//	//queue with visited bitset
//	u64 run_bfs_v5(const Idx& root){
//		setVPropAll(0xFFFF'FFFF'FFFF'FFFFUL);
//		const u64 numV = vArray.vSize;
//
//		Idx* __restrict activeList = (Idx*)globalAllocator.allocate(numV * sizeof(Idx));
//		vector<bool> visited(numV, 0);		// records visited vertices
//		u64 qHead = 1;
//		u64 qTail = 0;
//
//		activeList[0] = root;
//		vArray.setVProp(root, 0);
//
//		U64 iter = 1;
//		while(qHead - qTail) {
//			const U64 qHeadOrig = qHead;
//			#pragma omp parallel for num_threads(32)
//			for(U64 qT = qTail; qT < qHeadOrig; qT++){
//				//pop qTail
//				const u64 vIdx = activeList[qT];
//				const u64 outDeg = vArray.outDegree[vIdx];
//				const Idx* __restrict dstIdArr = vArray.dstIdArr[vIdx];
//				for(U64 i = 0; i < outDeg; i++){
//					const I64 dstId = dstIdArr[i];
//					if(!visited[dstId]){
//						visited[dstId] = 1;
//						vArray.prop[dstId] = iter;
//						U64 loc;
//						#pragma omp atomic capture
//						loc = qHead++;
//						activeList[loc] = dstId;
//					}
//				}
//			}
//			qTail = qHeadOrig;
//			iter++;
//		}
//		globalAllocator.deallocate(activeList, vArray.vSize * sizeof(Idx));
//		return iter;
//	}
//
//
//	// Same as stinger (incorrect version)
//	u64 run_page_rank_v1(VProp epsilon, VProp dampingfactor, u64 maxIter){		//using outgoing edges (push)
//
//		setVPropAll(1.0 / vArray.vSize);
//		const VProp offset = (1.0 - dampingfactor) / vArray.vSize;
//
//		VProp* vTemp = (VProp*)globalAllocator.allocate(vArray.vSize * sizeof(VProp));
//
//		u64 iter = 0;
//		VProp delta = 1.0;
//
//		while(delta > epsilon && iter != maxIter){
//
//			delta = 0.0;
//
//			//processing phase
//			#pragma omp parallel for reduction(+:delta)
//			for(u64 v = 0; v < vArray.vSize; v++){
//				vTemp[v] = 0;
//				const u64 outDeg = vArray.outDegree[v];
//				const Idx* __restrict dstIdArr = vArray.dstIdArr[v];
//				for(u64 e = 0; e < outDeg; e++){
//					const Idx dstId = dstIdArr[e];
//					//if(vArray.outDegree[dstId]){	//outDegree of dst can be zero
//						vTemp[v] += (vArray.prop[dstId] / vArray.outDegree[dstId]);
//					//}
//				}
//				vTemp[v] = vTemp[v] * dampingfactor + offset;
//				VProp myDelta = vTemp[v] - vArray.prop[v];
//				if(myDelta < 0){
//					myDelta = -myDelta;
//				}
//				delta += myDelta;
//			}
//
//			//apply phase
//			#pragma omp parallel for
//			for(u64 v = 0; v < vArray.vSize; v++){
//				vArray.prop[v] = vTemp[v];
//			}
//
//			iter++;
//		}
//
//		return iter;
//	}
//
//
//	// correct version
//	u64 run_page_rank_v2(VProp epsilon, VProp dampingfactor, u64 maxIter){		//using outgoing edges (push)
//
//		setVPropAll(1.0 / vArray.vSize);
//
//		VProp* vTemp = (VProp*)globalAllocator.allocate(vArray.vSize * sizeof(VProp));
//		#pragma omp parallel for
//		for(u64 v = 0; v < vArray.vSize; v++){
//			vTemp[v] = 0;
//		}
//
//
//		u64 iter = 0;
//
//		while(iter != maxIter){
//
//			//processing phase
//			#pragma omp parallel for
//			for(u64 v = 0; v < vArray.vSize; v++){
//				const u64 outDeg = vArray.outDegree[v];
//				if(outDeg){
//					const VProp pushVal = vArray.prop[v] / outDeg;
//					const Idx* __restrict dstIdArr = vArray.dstIdArr[v];
//					for(u64 e = 0; e < outDeg; e++){
//						const u64 dstId = dstIdArr[e];
//						#pragma omp atomic
//						vTemp[dstId] += pushVal;
//					}
//				}
//			}
//
//			//apply phase
//			#pragma omp parallel for
//			for(u64 v = 0; v < vArray.vSize; v++){
//				vArray.prop[v] = dampingfactor + (1 - dampingfactor) * vTemp[v];
//			}
//
//			iter++;
//		}
//
//		return iter;
//	}
//
//
//	u64 run_connected_components_v1(){
//		const u64 numV = vArray.vSize;
//		VProp* __restrict vProps = vArray.prop;
//
//		#pragma omp parallel for
//		for(u64 v = 0; v < numV; v++){
//			vProps[v] = v;
//		}
//
//		while(1){
//			bool isChanged = false;
//			#pragma omp parallel for
//			for(u64 v = 0; v < numV; v++){
//				const u64 numE = vArray.outDegree[v];
//				const Idx* __restrict dstIdArr = vArray.dstIdArr[v];
//				for(u64 e = 0; e < numE; e++){
//					const Idx dstId = dstIdArr[e];
//					if(vProps[v] > vProps[dstId]){
//						vProps[v] = vProps[dstId];
//						isChanged = true;
//					}
//				}
//			}
//
//			if(!isChanged){
//				break;
//			}
//
//			#pragma omp parallel for
//			for(u64 v = 0; v < numV; v++){
//				while(vProps[v] != vProps[vProps[v]]){
//					vProps[v] = vProps[vProps[v]];
//				}
//			}
//		}
//
//		u64 components = 1;
//		#pragma omp parallel for reduction(+:components)
//		for(u64 v = 1; v < numV; v++){
//			if(vProps[v] == v){
//				components++;
//			}
//		}
//
//		return components;
//	}
//
//
//#ifdef USE_INCOMING_EDGES
//	// correct version
//	u64 run_page_rank_v3(VProp epsilon, VProp dampingfactor, u64 maxIter){		//using incoming edges (pull)
//
//		setVPropAll(1.0 / vArray.vSize);
//		VProp* vTemp = (VProp*)globalAllocator.allocate(vArray.vSize * sizeof(VProp));
//
//		u64 iter = 0;
//
//		while(iter != maxIter){
//
//			#pragma omp parallel for
//			for(u64 v = 0; v < vArray.vSize; v++){
//				if(vArray.outDegree[v]){
//					vArray.prop[v] = vArray.prop[v] / vArray.outDegree[v];
//				}
//			}
//
//			#pragma omp parallel for
//			for(u64 v = 0; v < vArray.vSize; v++){
//				VProp sum = 0.0;
//				const u64 inDeg = vArray.inEdges[v].size();
//				const Idx* __restrict srcIdArr = vArray.inEdges[v].data();
//				for(u64 e = 0; e < inDeg; e++){
//					const Idx srcId = srcIdArr[e];
//					sum += vArray.prop[srcId];
//				}
//				vTemp[v] = dampingfactor + (1 - dampingfactor) * sum;
//			}
//
//			#pragma omp parallel for
//			for(u64 v = 0; v < vArray.vSize; v++){
//				vArray.prop[v] = vTemp[v];
//			}
//
//			iter++;
//		}
//
//		return iter;
//	}
//
//#endif
//
//
//
//	void setVPropAll(const VProp& vProp) {
//		#pragma omp parallel for
//		for(u64 v = 0; v < vArray.vSize; v++){
//			vArray.prop[v] = vProp;
//		}
//	}

};

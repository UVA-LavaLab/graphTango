#pragma once

#include <cstdint>
#include <immintrin.h>
#include <cassert>
#include <cstring>
#include <vector>
#include <map>

#include "common.h"
#include "CustomAllocator.h"

//typedef uint64_t u64;
//typedef uint32_t u32;
//typedef uint8_t  u8;

#define USE_64_BIT_KEY

#ifdef USE_64_BIT_KEY
typedef u64 Key;
#else
typedef u32 Key;
#endif

#define CACHE_LINE_SIZE			64U
#define ELEMS_IN_LINE			(CACHE_LINE_SIZE / sizeof(Key))

#ifdef USE_64_BIT_KEY
#define FLAG_EMPTY_SLOT			0xFFFFFFFFFFFFFFFFULL
#define FLAG_TOMB_STONE			0xFFFFFFFFFFFFFFFEULL
#define CONST_FACTOR_A			11400714818402812347ULL
#define NUM_INITIAL_ELEMS		(5 + 8)
#else
#define FLAG_EMPTY_SLOT			0xFFFFFFFFU
#define FLAG_TOMB_STONE			0xFFFFFFFEU
#define CONST_FACTOR_A			2654435769U
#endif



using namespace std;

#ifdef USE_CAHCE_FRIENDLY_HASH

template <typename Neigh>
class GraphiteHash{
public:

#ifdef USE_64_BIT_KEY
	static constexpr u8 rotation[8] = {0, 5, 2, 6, 1, 4, 7, 3};
	constexpr u32 getShiftAmt(){
		return 64 - __builtin_ctzl(capacity / ELEMS_IN_LINE);
	}
#else
	//static constexpr u8 rotation[16] = {0, 1, 4, 9, 15, 2, 10, 6, 3, 8, 14, 5, 12, 7, 11, 13};
	static constexpr u8 rotation[16] = {0, 9, 6, 15, 2, 11, 4, 13, 3, 7, 10, 14, 5, 1, 8, 12};
	constexpr u32 getShiftAmt(){
		return 32 - __builtin_ctzl(capacity / ELEMS_IN_LINE);
	}
#endif

	u32 degree = 0;
	u32 capacity = 0;
	Neigh* __restrict neighArr = nullptr;
	Neigh neigh[NUM_INITIAL_ELEMS];

	//map<u32, u32>	probes;

	GraphiteHash(){
		//capacity = NUM_INITIAL_ELEMS;
		//neighArr = (Neigh*)globalAllocator.allocate(capacity * sizeof(Neigh));
		//memset(neighArr, 0xff, capacity * sizeof(Neigh));
	}

	~GraphiteHash(){
		//free(neighArr);
	}

//	inline u32 hash1(u32 key){
//		return (u32)(key * CONST_FACTOR_A) >> (32 - capacity / ELEMS_IN_LINE + 1);
//	}
//
//	inline u32 hash2(u32 key){
//		const u32 mask = capacity / ELEMS_IN_LINE - 1;
//		return (key & mask) | 1;		//return an odd value
//	}
//
//	u32 find(u32 key){
//		const u32 numCacheLines = capacity / ELEMS_IN_LINE;
//		const u32 h1 = hash1(key);
//		for(u32 i = 0; i < numCacheLines; i++){
//			const u32 offset = ((h1 + i * hash2(key)) & (numCacheLines - 1)) * ELEMS_IN_LINE;
//			//cout << "\t" << offset << endl;
//
//			//check all of the brought cache line
//			//#pragma GCC unroll 8
//			for(int j = 0; j < ELEMS_IN_LINE; j++){
//				if(neighArr[offset + j] == key){
//					return offset + j;			//found
//				}
//				if(neighArr[offset + j] == FLAG_EMPTY_SLOT){
//					return FLAG_EMPTY_SLOT;		//not found
//				}
//			}
//		}
//		assert(false);		//should never reach here
//		return FLAG_EMPTY_SLOT;
//	}


	inline void insertDuringRehash(u32 key){
		const Key cacheLineMask = capacity / ELEMS_IN_LINE - 1;
		const Key h1 = (Key)(key * CONST_FACTOR_A) >> getShiftAmt();		//[0,1,...,#cache_lines]
		for(u32 i = 0; i <= cacheLineMask; i++){
			const Key h2 = (key & cacheLineMask) | 1;
			const Key base = ((h1 + i * h2) & cacheLineMask) * ELEMS_IN_LINE;	//cyclic within [0,1,...,#cache_lines]

			//check all elements of the cache line
			//#pragma GCC unroll 8
			for(int j = 0; j < ELEMS_IN_LINE; j++){
#ifdef USE_64_BIT_KEY
				const Key idx = base | ((key + rotation[j]) & 0x7);
#else
				const Key idx = base | ((key + rotation[j]) & 0xf);
#endif
				if(neighArr[idx].node == FLAG_EMPTY_SLOT){
					neighArr[idx].node = key;	//successful insertion
					return;
				}
			}
		}
	}

	void rehash(){
		Neigh* __restrict oldArr = neighArr;
		const u32 oldCap = capacity;

		capacity = capacity * 2;
		neighArr = (Neigh*)globalAllocator.allocPow2(capacity * sizeof(Neigh));
		memset(neighArr, 0xff, capacity * sizeof(Key));	//reset new array
		for(u32 i = 0; i < oldCap; i++){
			const Key key = oldArr[i].node;
			if(key < FLAG_TOMB_STONE){
				insertDuringRehash(key);
			}
		}

		globalAllocator.freePow2(oldArr, oldCap);
	}

	void insert(Key key, u64& edgeCnt){
		if(degree < HYBRID_HASH_PARTITION){
			//linear search
			for(u32 i = 0; i < degree; i++){
				if(neighArr[i].node == key){
					//found duplicate, nothing to do
					return;
				}
			}
			neighArr[degree].node = key;
			degree++;
			edgeCnt++;

			if(__builtin_expect(degree == capacity, 0)){
				// Two things can happen now.
				// 1. If we reached the partition threshold, switch to hash table
				// 2. Otherwise, just grow neighArr
				Neigh* __restrict oldArr = neighArr;
				const u32 oldCap = capacity;
				if(__builtin_expect(degree == HYBRID_HASH_PARTITION, 0)){
					//switch to hash table
					capacity = capacity * 4;
					neighArr = (Neigh*)globalAllocator.allocPow2(capacity * sizeof(Neigh));
					memset(neighArr, 0xff, capacity * sizeof(Neigh));	//reset new array
					for(u32 i = 0; i < degree; i++){
						insertDuringRehash(oldArr[i].node);
					}
				}
				else{
					//grow neighArr
					capacity = getNextPow2MinRet(capacity * 2);
					neighArr = (Neigh*)globalAllocator.allocPow2(capacity * sizeof(Neigh));
					memcpy(neighArr, oldArr, degree * sizeof(Neigh));
				}
				if(oldArr != neigh){
					globalAllocator.freePow2(oldArr, oldCap);
				}
			}
		}
		else{
			//hash based search
			if(__builtin_expect(degree > (capacity / 2), 0)){
				//Load factor is 0.5. Grow table
				rehash();
			}

			//u32 probelen = 0;
			const Key cacheLineMask = capacity / ELEMS_IN_LINE - 1;
			const Key h1 = (Key)(key * CONST_FACTOR_A) >> getShiftAmt();		//[0,1,...,#cache_lines]
			for(Key i = 0; i <= cacheLineMask; i++){
				const Key h2 = (key & cacheLineMask) | 1;
				const Key base = ((h1 + i * h2) & cacheLineMask) * ELEMS_IN_LINE;	//cyclic within [0,1,...,#cache_lines]

				//check all elements of the cache line
				//#pragma GCC unroll 8
				for(int j = 0; j < ELEMS_IN_LINE; j++){
	#ifdef USE_64_BIT_KEY
					const Key idx = base | ((key + rotation[j]) & 0x7);
	#else
					const Key idx = base | ((key + rotation[j]) & 0xf);
	#endif
					//probelen++;
					if(neighArr[idx].node >= FLAG_TOMB_STONE){
						neighArr[idx].node = key;	//successful insertion
						degree++;
						edgeCnt++;
						//probes[probelen]++;
						return;
					}
					if(neighArr[idx].node == key){
						return;					//found, no need to do anything
					}
				}
			}
		}
	}

	void erase(u32 key){

	}

};

template<typename Neigh>
constexpr u8 GraphiteHash<Neigh>::rotation[8];

#endif


#ifdef USE_CAHCE_FRIENDLY_HASH_ONLY

template <typename Neigh>
class GraphiteHash{
public:

#ifdef USE_64_BIT_KEY
	static constexpr u8 rotation[8] = {0, 5, 2, 6, 1, 4, 7, 3};
	constexpr u32 getShiftAmt(){
		return 64 - __builtin_ctzl(capacity / ELEMS_IN_LINE);
	}
#else
	//static constexpr u8 rotation[16] = {0, 1, 4, 9, 15, 2, 10, 6, 3, 8, 14, 5, 12, 7, 11, 13};
	static constexpr u8 rotation[16] = {0, 9, 6, 15, 2, 11, 4, 13, 3, 7, 10, 14, 5, 1, 8, 12};
	constexpr u32 getShiftAmt(){
		return 32 - __builtin_ctzl(capacity / ELEMS_IN_LINE);
	}
#endif

	u32 degree = 0;
	u32 capacity = 0;
	Neigh* __restrict neighArr = nullptr;
	Neigh* __restrict adjList = nullptr;
	Neigh neigh[NUM_INITIAL_ELEMS];

	//map<u32, u32>	probes;

	GraphiteHash(){
		//capacity = NUM_INITIAL_ELEMS;
		//neighArr = (Neigh*)globalAllocator.allocate(capacity * sizeof(Neigh));
		//memset(neighArr, 0xff, capacity * sizeof(Neigh));
	}

	~GraphiteHash(){
		//free(neighArr);
	}

	inline void insertDuringRehash(u32 key){
		const Key cacheLineMask = capacity / ELEMS_IN_LINE - 1;
		const Key h1 = (Key)(key * CONST_FACTOR_A) >> getShiftAmt();		//[0,1,...,#cache_lines]
		for(u32 i = 0; i <= cacheLineMask; i++){
			const Key h2 = (key & cacheLineMask) | 1;
			const Key base = ((h1 + i * h2) & cacheLineMask) * ELEMS_IN_LINE;	//cyclic within [0,1,...,#cache_lines]

			//check all elements of the cache line
			//#pragma GCC unroll 8
			for(int j = 0; j < ELEMS_IN_LINE; j++){
#ifdef USE_64_BIT_KEY
				const Key idx = base | ((key + rotation[j]) & 0x7);
#else
				const Key idx = base | ((key + rotation[j]) & 0xf);
#endif
				if(neighArr[idx].node == FLAG_EMPTY_SLOT){
					neighArr[idx].node = key;	//successful insertion
					return;
				}
			}
		}
	}

	void rehash(){
		Neigh* __restrict oldArr = neighArr;
		const u32 oldCap = capacity;

		capacity = capacity * 2;
		neighArr = (Neigh*)globalAllocator.allocPow2(capacity * sizeof(Neigh));
		memset(neighArr, 0xff, capacity * sizeof(Key));	//reset new array
		for(u32 i = 0; i < oldCap; i++){
			const Key key = oldArr[i].node;
			if(key < FLAG_TOMB_STONE){
				insertDuringRehash(key);
			}
		}

		globalAllocator.freePow2(oldArr, oldCap);
	}

	void insert(Key key, u64& edgeCnt){
		if(degree == NUM_INITIAL_ELEMS){
			//switch to hash table, nut key is not yet inserted
			capacity = getNextPow2MinRet(NUM_INITIAL_ELEMS * 4);
			neighArr = (Neigh*)globalAllocator.allocPow2(capacity * sizeof(Neigh));
			memset(neighArr, 0xff, capacity * sizeof(Neigh));	//reset new array
			for(u32 i = 0; i < NUM_INITIAL_ELEMS; i++){
				insertDuringRehash(neigh[i].node);
			}
		}

		if(degree < NUM_INITIAL_ELEMS){
			//linear search
			for(u32 i = 0; i < degree; i++){
				if(neigh[i].node == key){
					//found duplicate, nothing to do
					return;
				}
			}
			neigh[degree].node = key;
			degree++;
			edgeCnt++;
		}
		else{
			//hash based search
			if(__builtin_expect(degree > (capacity / 2), 0)){
				//Load factor is 0.5. Grow table
				rehash();
			}

			//u32 probelen = 0;
			const Key cacheLineMask = capacity / ELEMS_IN_LINE - 1;
			const Key h1 = (Key)(key * CONST_FACTOR_A) >> getShiftAmt();		//[0,1,...,#cache_lines]
			for(Key i = 0; i <= cacheLineMask; i++){
				const Key h2 = (key & cacheLineMask) | 1;
				const Key base = ((h1 + i * h2) & cacheLineMask) * ELEMS_IN_LINE;	//cyclic within [0,1,...,#cache_lines]

				//check all elements of the cache line
				//#pragma GCC unroll 8
				for(int j = 0; j < ELEMS_IN_LINE; j++){
	#ifdef USE_64_BIT_KEY
					const Key idx = base | ((key + rotation[j]) & 0x7);
	#else
					const Key idx = base | ((key + rotation[j]) & 0xf);
	#endif
					//probelen++;
					if(neighArr[idx].node >= FLAG_TOMB_STONE){
						neighArr[idx].node = key;	//successful insertion
						degree++;
						edgeCnt++;
						//probes[probelen]++;
						adjList = nullptr;	//adjacency list no longer valid
						return;
					}
					if(neighArr[idx].node == key){
						return;					//found, no need to do anything
					}
				}
			}
		}
	}

	void erase(Key key, u64& edgeCnt){

	}

};

template<typename Neigh>
constexpr u8 GraphiteHash<Neigh>::rotation[8];

#endif





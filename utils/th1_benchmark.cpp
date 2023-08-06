#include <cstdint>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <cstring>
#include <chrono>
#include <algorithm>
#include <random>

typedef uint8_t u8;
typedef uint32_t u32;
typedef uint64_t u64;

using namespace std;
using namespace std::chrono;

#define CACHE_LINE_SIZE				64
#define NUM_TOT_ARRAY_ELEMS			(1024ULL * 1024ULL * 128)
#define TOT_ARRAY_SIZE				(CACHE_LINE_SIZE * NUM_TOT_ARRAY_ELEMS + 1024*CACHE_LINE_SIZE)
#define U32_ELEM_PER_CACHE_LINE		(CACHE_LINE_SIZE/sizeof(u32))

#define NUM_TRIALS					40000000ULL

class CacheLine {
public:
	u32 pad[U32_ELEM_PER_CACHE_LINE];
	
	//CacheLine(){
	//	for(u64 i = 0; i < U32_ELEM_PER_CACHE_LINE; i++){
	//		pad[i] = rand();
	//	}
	//}
};

int main(int argc, const char** argv){
	if(argc != 2){
		cerr << "usage: " << argv[0] << " <edge_elem_size>" << endl;
		exit(-1);
	}

	srand(42);
	u64 elemSize = atol(argv[1]);
	u64 elemPerCacheLine = CACHE_LINE_SIZE / elemSize;

	assert(sizeof(CacheLine) == CACHE_LINE_SIZE);
	CacheLine* __restrict arr = (CacheLine*)aligned_alloc(CACHE_LINE_SIZE, TOT_ARRAY_SIZE);
	assert(arr);
	memset(arr, 0, TOT_ARRAY_SIZE);
	u32* __restrict idxArr = (u32*)aligned_alloc(CACHE_LINE_SIZE, NUM_TOT_ARRAY_ELEMS*sizeof(u32));
	u32* __restrict idxPtr = idxArr;
	assert(idxArr);
	for(u64 i = 0; i < NUM_TOT_ARRAY_ELEMS; i++){
		idxArr[i] = i;
	}
	shuffle(idxArr, idxArr + NUM_TOT_ARRAY_ELEMS, default_random_engine(47));
	u64 idxArrIdx = 0;

	u32 garbage = 0;
	double hashTime = 0.0;

	// hash map access simulation
	{
		idxArr = idxPtr;
		auto start = high_resolution_clock::now();
		for(u64 i = 0; i < NUM_TRIALS; i++){
			u32 idx = *idxArr++;
			garbage += arr[idx].pad[0];
			
			idx = *idxArr++;
			garbage += arr[idx].pad[1];
			
			idx = *idxArr++;
			garbage += arr[idx].pad[2];
		}
		auto end = high_resolution_clock::now();
		hashTime = (1e-6)*(end-start).count();
		cout << "Hashmap time (ms): " << hashTime << endl;
		//cout << "garbage: " << garbage << endl;
	}

	u64 lineOffset = 0;
	while(true){
		idxArr = idxPtr;
		u64 extraLine = 1ULL << lineOffset;
		auto start = high_resolution_clock::now();
		for(u64 i = 0; i < NUM_TRIALS; i++){
			u32 idx = *idxArr++;
			garbage += arr[idx].pad[0];
			
			idx = *idxArr++;
			for(u64 j = 0; j < extraLine; j++){
				garbage += arr[idx+j].pad[0];
			}
		}
		
		auto end = high_resolution_clock::now();
		double adjTime = (1e-6)*(end-start).count();
		cout << "Adjacency list time [TH1 = " << extraLine * elemPerCacheLine << "] " << adjTime << endl;
		if(adjTime > hashTime) break;
		lineOffset++;
	}
	lineOffset--;
	
	cout << "Recommended TH1: " << (1ULL<<lineOffset) * elemPerCacheLine << endl;
	cout << garbage << endl;
	
	return 0;
}

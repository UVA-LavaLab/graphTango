#pragma once

#include <omp.h>
#include <cstdlib>
#include <iostream>
#include <vector>
#include "common.h"

#define	BLOCK_BITS			20
#define BLOCK_SIZE			(1UL << BLOCK_BITS)

template <u64 MAX_THREADS = 32, u64 MAX_SEGMENT_BITS = 32>
class LockFreePool {
	alignas(64) std::vector<void* __restrict> freeSegments[MAX_THREADS][MAX_SEGMENT_BITS];

public:

	void* allocPow2(u64 size){
		//return segments[omp_get_thread_num()][getPow2Log2(size)].alloc(size);
		auto& vec = freeSegments[omp_get_thread_num()][getPow2Log2(size)];
		if(__builtin_expect(vec.size() == 0, 0)){
			//create and append a block
			const u64 blockSize = (size <= BLOCK_SIZE) ? BLOCK_SIZE : size;
			u8* __restrict const  blockStart = (u8*)aligned_alloc(64, blockSize);
			u8* __restrict blockEnd = blockStart + blockSize;
			//blocks.push_back(blockStart);
			const u64 numNewSegments = blockSize >> getPow2Log2(size);
			vec.reserve(vec.size() + numNewSegments);
			while(blockStart != blockEnd){
				vec.push_back(blockEnd);
				blockEnd -= size;
			}
		}
		void* ret = vec.back();
		vec.pop_back();
		return ret;
	}

	void* alloc(u64 size){
		if(size == 0){
			return nullptr;
		}
		return allocPow2(getNextPow2(size));
	}

	void freePow2(void* ptr, u64 size){
		freeSegments[omp_get_thread_num()][getPow2Log2(size)].push_back(ptr);
	}

	void free(void* ptr, u64 size){
		freePow2(ptr, getNextPow2(size));
	}
};


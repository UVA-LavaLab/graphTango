#pragma once

#include <cstdint>
#include <cstdio>
#include "GEmpty.h"

typedef 	uint64_t	U64;
typedef		int64_t		I64;
typedef		uint32_t	U32;
typedef		int64_t		I32;
typedef		uint8_t		U8;
typedef		int8_t		I8;

typedef 	uint64_t	u64;
typedef		int64_t		i64;
typedef		uint32_t	u32;
typedef		int64_t		i32;
typedef		uint8_t		u8;
typedef		int8_t		i8;

typedef		I64			Idx;

#define 	MIN_RET_VAL  		1

//define only one of the following
//#define 	USE_HYBRID_HASHMAP
#define 	USE_HYBRID_HASHMAP_WITH_CFH
//#define		USE_ONLY_LINEAR
//#define		USE_ONLY_HASHMAP
//#define		USE_WEIRD_SCHEME
//#define		USE_HYBRID_HASHMAP_WITH_GROUPING
//#define		USE_HYBRID_HASHMAP_WITH_GROUPING_TIGHTER
//#define		USE_HYBRID_HASHMAP_WITH_GROUPING_AND_EDGE_ARR_LOCKING
//#define		USE_SORTED_EDGES
//#define		USE_CAHCE_FRIENDLY_HASH
//#define		USE_CAHCE_FRIENDLY_HASH_ONLY

#if			defined(USE_HYBRID_HASHMAP)													\
			|| defined(USE_HYBRID_HASHMAP_WITH_CFH)										\
			|| defined(USE_HYBRID_HASHMAP_WITH_GROUPING)								\
			|| defined(USE_HYBRID_HASHMAP_WITH_GROUPING_AND_EDGE_ARR_LOCKING)			\
			|| defined(USE_HYBRID_HASHMAP_WITH_GROUPING_TIGHTER)						\
			|| defined(USE_CAHCE_FRIENDLY_HASH)
#define		HYBRID_HASH_PARTITION		64UL
#endif

#ifdef		USE_SORTED_EDGES
#define		LINEAR_BUFF_SIZE			512UL
#endif

typedef struct {
	u32 dst;
	u32 loc;
} DstLocPair;

// Log2 for power of 2 integers
//#define 	LOG2(x) 	__builtin_ctzl(x)
constexpr U64 getPow2Log2(U64 val) {
	return __builtin_ctzl(val);
}

static U64 getNextPow2(U64 val) {
	if(val == 0){
		return 1;
	}
	U64 res = 1UL << (63 - __builtin_clzl(val));
	if(res < val){
		res = res << 1;
	}
	return res;
}

static U64 getNextPow2MinRet(U64 val) {
	if(val <= MIN_RET_VAL){
		return MIN_RET_VAL;
	}
	U64 res = 1UL << (63 - __builtin_clzl(val));
	if(res < val){
		res = res << 1;
	}
	return res;
}


constexpr U64 getNextPow2Log2(U64 val) {
	//return 1 << (64 - __builtin_clzl(val - 1));
	return (64 - __builtin_clzl(val - 1));
}

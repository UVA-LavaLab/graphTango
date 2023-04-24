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
typedef		int32_t		i32;
typedef		uint8_t		u8;
typedef		int8_t		i8;

typedef		I64			Idx;

#define 	MIN_RET_VAL  		2

//#define LIKWID_PERFMON

// This block enables to compile the code with and without the likwid header in place
#ifdef LIKWID_PERFMON
#include <likwid-marker.h>
#else
#define LIKWID_MARKER_INIT
#define LIKWID_MARKER_THREADINIT
#define LIKWID_MARKER_SWITCH
#define LIKWID_MARKER_REGISTER(regionTag)
#define LIKWID_MARKER_START(regionTag)
#define LIKWID_MARKER_STOP(regionTag)
#define LIKWID_MARKER_CLOSE
#define LIKWID_MARKER_GET(regionTag, nevents, events, time, count)
#endif

//#define		CALC_TYPE_SWITCH
//#define		USE_CFH_FOR_DAH
//#define		CALC_MEM_PER_EDGE
//#define		ALLOC_FROM_NEXT_BLOCK
//#define 	CALC_EDGE_TOUCHED
//#define 	USE_HUGEPAGE

//define only one of the following
//#define 	USE_HYBRID_HASHMAP
//#define 	USE_HYBRID_HASHMAP_WITH_CFH
#define 	USE_GT_LOAD_BALANCED
//#define 	USE_GT_BALANCED
//#define 	USE_GT_BALANCED_TYPE3_ONLY
//#define 	USE_GT_BALANCED_MALLOC
//#define 	USE_GT_BALANCED_STDMAP
//#define 	USE_GT_BALANCED_ABSEIL
//#define 	USE_GT_BALANCED_RHH
//#define 	USE_GT_BALANCED_TSL_RHH
//#define 	USE_GT_BALANCED_MALLOC_STDMAP
//#define 	USE_GT_BALANCED_DYN_PARTITION
//#define 	USE_GT_UPDATE
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
			|| defined(USE_CAHCE_FRIENDLY_HASH)											\
			|| defined(USE_GT_BALANCED)													\
			|| defined(USE_GT_UPDATE)													\
			|| defined(USE_GT_BALANCED_MALLOC)											\
			|| defined(USE_GT_BALANCED_STDMAP)											\
			|| defined(USE_GT_BALANCED_MALLOC_STDMAP)									\
			|| defined(USE_GT_BALANCED_DYN_PARTITION)									\
			|| defined(USE_GT_BALANCED_ABSEIL)											\
			|| defined(USE_GT_BALANCED_RHH)												\
			|| defined(USE_GT_BALANCED_TSL_RHH)											\
			|| defined(USE_GT_LOAD_BALANCED)
#define		HYBRID_HASH_PARTITION		64UL
#endif

#ifdef		USE_SORTED_EDGES
#define		LINEAR_BUFF_SIZE			512UL
#endif

#ifdef USE_GT_LOAD_BALANCED
#define		LB_NUMBER_OF_BUCKETS		256UL
#endif

typedef struct {
	u32 dst;
	u32 loc;
} DstLocPair;

typedef enum {
	VTYPE_1,
	VTYPE_2,
	VTYPE_3
} VType;

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

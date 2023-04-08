
//#include "MemoryPool.h"

//MemoryPool globalAllocator;

#include "LockFreePoolWithList.h"

#if defined(USE_GT_BALANCED_MALLOC) || defined(USE_GT_BALANCED_MALLOC_STDMAP)
#else
LockFreePoolWithList<> globalAllocator;
#endif

#ifdef CALC_EDGE_TOUCHED
u64 g_edge_touched = 0;
#endif

#ifdef CALC_PROBING_DISTANCE
u64 g_num_type3 = 0;
u64 g_num_probes = 0;
#endif


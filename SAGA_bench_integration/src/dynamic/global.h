#pragma once

#include "LockFreePoolWithList.h"

#if defined(USE_GT_BALANCED_MALLOC) || defined(USE_GT_BALANCED_MALLOC_STDMAP)
#else
extern LockFreePoolWithList<> globalAllocator;
#endif

#ifdef CALC_EDGE_TOUCHED
extern u64 g_edge_touched;
#endif



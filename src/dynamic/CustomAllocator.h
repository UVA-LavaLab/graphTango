#pragma once

#if defined(USE_GT_BALANCED_MALLOC) || defined(USE_GT_BALANCED_MALLOC_STDMAP)
#else

#include "LockFreePoolWithList.h"
#include "global.h"

#include <cstdlib>
#include <memory>
#include <cassert>

template <class T>
struct custom_allocator {
  typedef T value_type;

  custom_allocator() noexcept {}

  template <class U> custom_allocator (const custom_allocator<U>&) noexcept {}

  T* allocate (std::size_t n) {
	  return static_cast<T*>(globalAllocator.allocate(n * sizeof(T)));
  }

  void deallocate (T* p, std::size_t n) {
	  globalAllocator.deallocate(p, n * sizeof(T));
  }
};

template <class T, class U>
constexpr bool operator== (const custom_allocator<T>&, const custom_allocator<U>&) noexcept
{return true;}

template <class T, class U>
constexpr bool operator!= (const custom_allocator<T>&, const custom_allocator<U>&) noexcept
{return false;}

#endif

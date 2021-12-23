#ifndef SLIDING_QUEUE_H_
#define SLIDING_QUEUE_H_

#include <algorithm>
#include <cassert>

#include "../common/platform_atomics.h"

/* This implementation has been borrowed from GAP Benchmark Suite (https://github.com/sbeamer/gapbs) 
   and modified for SAGA-Bench. */

template <typename T>
class QueueBuffer;

template <typename T>
class SlidingQueue {
  T *shared;
  size_t shared_in;
  size_t shared_out_start;
  size_t shared_out_end;
  T* curr;
  T* next;
  // size_t _shared_size;
  friend class QueueBuffer<T>;

 public:
  explicit SlidingQueue(size_t shared_size) {    
    curr = new T[shared_size];
    next = new T[shared_size];
    shared = curr;
    reset();
  }

  ~SlidingQueue() {
    delete[] curr;
    delete[] next;   
  }

  void push_back(T to_add) {
    next[shared_in++] = to_add;    
  }

  bool empty() const {
    return shared_out_start == shared_out_end;
  }

  void reset() {
    shared_out_start = 0;
    shared_out_end = 0;
    shared_in = 0;
  }

  void slide_window() {    
    shared_out_start = 0;
    shared_out_end = shared_in;
    shared_in = 0;
    shared = next;
    next = curr;
    curr = shared;
  }

  typedef T* iterator;

  iterator begin() const {
    return shared + shared_out_start;
  }

  iterator end() const {
    return shared + shared_out_end;
  }

  size_t size() const {
    return end() - begin();
  }
};


template <typename T>
class QueueBuffer {
  size_t in;
  T *local_queue;
  SlidingQueue<T> &sq;
  const size_t local_size;

 public:
  explicit QueueBuffer(SlidingQueue<T> &master, size_t given_size = 16384)
      : sq(master), local_size(given_size) {
    in = 0;
    local_queue = new T[local_size];
  }

  ~QueueBuffer() {
    delete[] local_queue;
  }

  void push_back(T to_add) {
    if (in == local_size)
      flush();
    local_queue[in++] = to_add;
  }

  void flush() {
    T* shared_queue = sq.next;
    //T *shared_queue = sq.shared;
    size_t copy_start = fetch_and_add(sq.shared_in, in);
    std::copy(local_queue, local_queue+in, shared_queue+copy_start);
    in = 0;
  }
};

#endif  // SLIDING_QUEUE_H_

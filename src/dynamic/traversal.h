#ifndef TRAVERSAL_H_
#define TRAVERSAL_H_

#include "types.h"
#include "adListShared.h"
#include "stinger.h"
#include "darhh.h"
#include "adListChunked.h"
#include "GraphTango.h"
#include "Vertex.h"

#include "topDataStruc.h"

template<typename T>
class neighborhood;

template<typename T>
class neighborhood_iter {
public:
	neighborhood_iter(T *ds, NodeID n, bool in_neigh);
	bool operator!=(neighborhood_iter &it);
	neighborhood_iter& operator++();
	neighborhood_iter& operator++(int);
	NodeID operator*();
	Weight extractWeight();
};

template<typename U>
class neighborhood_iter<GraphTango<U>> {
	friend class neighborhood<GraphTango<U>> ;
private:
	u32 idx;
	U* arr;
	U** blocks;
	bool isType3;
public:
	neighborhood_iter(u32 idx, U* arr, U** blocks, bool isType3) :
			idx(idx),
			arr(arr),
			blocks(blocks),
			isType3(isType3) { }

	bool operator!=(const neighborhood_iter<GraphTango<U>> &it) {
		return idx != it.idx;
	}

	neighborhood_iter& operator++() {
#ifdef CALC_EDGE_TOUCHED
		#pragma omp atomic
		g_edge_touched++;
#endif
		idx++;
		if(isType3 && !(idx % EdgeArray<U>::BLOCK_SIZE)){
			blocks++;
			arr = *blocks;
		}
		return *this;
	}

	neighborhood_iter& operator++(int) {
#ifdef CALC_EDGE_TOUCHED
		#pragma omp atomic
		g_edge_touched++;
#endif
		idx++;
		if(isType3 && !(idx % EdgeArray<U>::BLOCK_SIZE)){
			blocks++;
			arr = *blocks;
		}
		return *this;
	}

	NodeID operator*() const {
		u32 offset = idx % EdgeArray<U>::BLOCK_SIZE;
		return arr[offset].getNodeID();
	}

	Weight extractWeight() const {
		u32 offset = idx % EdgeArray<U>::BLOCK_SIZE;
		return arr[offset].getWeight();
	}
};

template<typename U>
class neighborhood_iter<adList<U>> {
	friend class neighborhood<adList<U>> ;
private:
	adList<U> *ds;
	NodeID node;
	bool in_neigh;
	U *cursor;
public:
	neighborhood_iter(adList<U> *_ds, NodeID _n, bool _in_neigh) :
			ds(_ds), node(_n), in_neigh(_in_neigh) {
		if (in_neigh) {
			bool empty = ds->in_neighbors[node].empty();
			cursor = empty ? 0 : &(ds->in_neighbors[node][0]);
		} else {
			bool empty = ds->out_neighbors[node].empty();
			cursor = empty ? 0 : &(ds->out_neighbors[node][0]);
		}
	}

	bool operator!=(const neighborhood_iter<adList<U>> &it) {
		return cursor != it.cursor;
	}

	neighborhood_iter& operator++() {
		if (in_neigh) {
			int size_in_neigh = ds->in_neighbors[node].size();
			if (cursor == &(ds->in_neighbors[node][size_in_neigh - 1])) {
				cursor = nullptr;
			} else
				cursor = cursor + 1;
		} else {
			int size_out_neigh = ds->out_neighbors[node].size();
			if (cursor == &(ds->out_neighbors[node][size_out_neigh - 1])) {
				cursor = nullptr;
			} else
				cursor = cursor + 1;
		}

		return *this;
	}

	neighborhood_iter& operator++(int) {
		if (in_neigh) {
			int size_in_neigh = ds->in_neighbors[node].size();
			if (cursor == &(ds->in_neighbors[node][size_in_neigh - 1])) {
				cursor = nullptr;
			} else
				cursor = cursor + 1;
		} else {
			int size_out_neigh = ds->out_neighbors[node].size();
			if (cursor == &(ds->out_neighbors[node][size_out_neigh - 1])) {
				cursor = nullptr;
			} else
				cursor = cursor + 1;
		}

		return *this;
	}

	NodeID operator*() {
		return cursor->getNodeID();
	}

	Weight extractWeight() {
		return cursor->getWeight();
	}
};

template<typename U>
class neighborhood_iter<adListShared<U>> {
	friend class neighborhood<adListShared<U>> ;
private:
	U *cursor;
	U *endPtr;
public:
	neighborhood_iter(adListShared<U> *ds, NodeID node, bool in_neigh) {
		cursor = nullptr;
		endPtr = nullptr;

		if(in_neigh){
			const int sz = ds->in_neighbors[node].size();
			if(sz){
				cursor = &(ds->in_neighbors[node][0]);
				endPtr = cursor + sz;
			}
		}
		else{
			const int sz = ds->out_neighbors[node].size();
			if(sz){
				cursor = &(ds->out_neighbors[node][0]);
				endPtr = cursor + sz;
			}
		}

		while(cursor < endPtr){
			if(cursor->node >= 0){
				//found a valid node
				return;
			}
			cursor++;
		}
		cursor = nullptr;
		endPtr = nullptr;
	}

	bool operator!=(const neighborhood_iter<adListShared<U>> &it) {
		return cursor != it.cursor;
	}

	neighborhood_iter& operator++() {
		cursor++;
		while(cursor < endPtr){
			if(cursor->node >= 0){
				//found a valid node
				return *this;
			}
			cursor++;
		}
		cursor = nullptr;
		return *this;
	}

	neighborhood_iter& operator++(int) {
		cursor++;
		while(cursor < endPtr){
			if(cursor->node >= 0){
				//found a valid node
				return *this;
			}
			cursor++;
		}
		cursor = nullptr;
		return *this;
	}

	NodeID operator*() {
		return cursor->getNodeID();
	}

	Weight extractWeight() {
		return cursor->getWeight();
	}
};

/*
//specialization for stinger ---- OLD
template<>
class neighborhood_iter<stinger> {
	friend class neighborhood<stinger> ;

private:
	stinger *ds;
	NodeID node;
	bool in_neigh;
	stinger_edge *cursor;
	stinger_eb *curr_eb;
	stinger_vertex *sv;
	int cursor_index;

public:
	neighborhood_iter(stinger *_ds, NodeID _n, bool _in_neigh) :
			ds(_ds), node(_n), in_neigh(_in_neigh) {
		sv = &(ds->vertices[node]);

		if (in_neigh) {
			bool empty = (sv->in_neighbors->numEdges == 0);
			cursor = empty ? 0 : &(sv->in_neighbors->edges[0]);
			curr_eb = sv->in_neighbors;
			if (!empty)
				cursor_index = 0;
		}

		else {
			bool empty = (sv->out_neighbors->numEdges == 0);
			cursor = empty ? 0 : &(sv->out_neighbors->edges[0]);
			curr_eb = sv->out_neighbors;
			if (!empty)
				cursor_index = 0;
		}
	}

	bool operator!=(const neighborhood_iter<stinger> &it) {
		return cursor != it.cursor;
	}

	neighborhood_iter& operator++() {
		// just increment by 1 if we are in an edgeblock and more edges left
		if (cursor_index < (curr_eb->numEdges - 1)) {
			cursor = cursor + 1;
			cursor_index++;
		}

		else if (cursor_index == (curr_eb->numEdges - 1)) {
			// We are done with current edgeblock
			if (curr_eb->next != nullptr) {
				// move to next one, if there is one
				curr_eb = curr_eb->next;
				cursor = &(curr_eb->edges[0]);
				cursor_index = 0;
			} else {
				// there is no further, end of traversal
				cursor = nullptr;
			}
		}

		return *this;
	}

	neighborhood_iter& operator++(int) {
		// just increment by 1 if we are in an edgeblock and more edges left
		if (cursor_index < (curr_eb->numEdges - 1)) {
			cursor = cursor + 1;
			cursor_index++;
		}

		else if (cursor_index == (curr_eb->numEdges - 1)) {
			// We are done with current edgeblock
			if (curr_eb->next != nullptr) {
				// move to next one, if there is one
				curr_eb = curr_eb->next;
				cursor = &(curr_eb->edges[0]);
				cursor_index = 0;
			} else {
				// there is no further, end of traversal
				cursor = nullptr;
			}
		}

		return *this;
	}

	NodeID operator*() {
		assert(cursor->neighbor != -1);
		return cursor->neighbor;
	}

	Weight extractWeight() {
		return cursor->weight;
	}
};
*/


//specialization for stinger ---- NEW
template<>
class neighborhood_iter<stinger> {
	friend class neighborhood<stinger> ;

private:
	//stinger *ds;
	//NodeID node;
	//bool in_neigh;
	stinger_edge *cursor;
	stinger_eb *curr_eb;
	int cursor_index;
	//stinger_vertex *sv;


public:
	neighborhood_iter(stinger *_ds, NodeID _n, bool _in_neigh) {
		stinger_vertex* sv = &(_ds->vertices[_n]);

		if (_in_neigh) {
			curr_eb = sv->in_neighbors;
		}
		else{
			curr_eb = sv->out_neighbors;
		}

		cursor_index = 0;

		//find first edge
		while(curr_eb){
			while(cursor_index < curr_eb->high){
				if(curr_eb->edges[cursor_index].neighbor >= 0){
					cursor = curr_eb->edges + cursor_index;
					return;
				}
				cursor_index++;
			}

			cursor_index = 0;
			curr_eb = curr_eb->next;
		}
		cursor = nullptr;
	}

	bool operator!=(const neighborhood_iter<stinger> &it) {
		return cursor != it.cursor;
	}

	neighborhood_iter& operator++() {
		while(curr_eb){
			cursor_index++;
			while(cursor_index < curr_eb->high){
				if(curr_eb->edges[cursor_index].neighbor >= 0){
					//valid edge (otherwise deleted)
					cursor = curr_eb->edges + cursor_index;
					return *this;
				}
				cursor_index++;
			}

			cursor_index = 0;
			curr_eb = curr_eb->next;
		}

		cursor = nullptr;
		return *this;
	}

	neighborhood_iter& operator++(int) {
		while(curr_eb){
			cursor_index++;
			while(cursor_index < curr_eb->high){
				if(curr_eb->edges[cursor_index].neighbor >= 0){
					//valid edge (otherwise deleted)
					cursor = curr_eb->edges + cursor_index;
					return *this;
				}
				cursor_index++;
			}

			cursor_index = 0;
			curr_eb = curr_eb->next;
		}

		cursor = nullptr;
		return *this;
	}

	NodeID operator*() {
		//assert(cursor->neighbor != -1);
		return cursor->neighbor;
	}

	Weight extractWeight() {
		return cursor->weight;
	}
};



//// // ------------------------------adList_chunk---------------------------------- OLD
//template<typename U>
//class neighborhood_iter<adListChunked<U>> {
//	friend class neighborhood<adListChunked<U>> ;
//private:
//	adListChunked<U> *ds;
//	NodeID node;
//	bool in_neigh;
//	U *cursor;
//
//	int64_t part_idx;
//	int64_t sub_idx;
//
//public:
//	neighborhood_iter(adListChunked<U> *_ds, NodeID _n, bool _in_neigh) :
//			ds(_ds), node(_n), in_neigh(_in_neigh) {
//		// part_idx = _n % (ds->num_partitions);
//		// sub_idx = (int) _n/(ds->num_partitions);
//		part_idx = ds->pt_hash(_n);
//		sub_idx = ds->hash_within_chunk(_n);
//		if (in_neigh) {
//			bool empty = ds->in[part_idx]->partAdList->neighbors[sub_idx].empty();
//			cursor = empty ? nullptr : &(ds->in[part_idx]->partAdList->neighbors[sub_idx][0]);
//		} else {
//			bool empty = ds->out[part_idx]->partAdList->neighbors[sub_idx].empty();
//			cursor = empty ? nullptr : &(ds->out[part_idx]->partAdList->neighbors[sub_idx][0]);
//		}
//	}
//
//	bool operator!=(const neighborhood_iter<adListChunked<U>> &it) {
//		return cursor != it.cursor;
//	}
//
//	neighborhood_iter& operator++() {
//		if (in_neigh) {
//			int size_in_neigh = ds->in[part_idx]->partAdList->neighbors[sub_idx].size();
//			if (cursor == &(ds->in[part_idx]->partAdList->neighbors[sub_idx][size_in_neigh - 1]))
//				cursor = nullptr;
//			else
//				cursor = cursor + 1;
//		} else {
//			int size_out_neigh = ds->out[part_idx]->partAdList->neighbors[sub_idx].size();
//			if (cursor == &(ds->out[part_idx]->partAdList->neighbors[sub_idx][size_out_neigh - 1]))
//				cursor = nullptr;
//			else
//				cursor = cursor + 1;
//		}
//
//		return *this;
//	}
//
//	neighborhood_iter& operator++(int) {
//		if (in_neigh) {
//			int size_in_neigh = ds->in[part_idx]->partAdList->neighbors[sub_idx].size();
//			if (cursor == &(ds->in[part_idx]->partAdList->neighbors[sub_idx][size_in_neigh - 1]))
//				cursor = nullptr;
//			else
//				cursor = cursor + 1;
//		} else {
//			int size_out_neigh = ds->out[part_idx]->partAdList->neighbors[sub_idx].size();
//			if (cursor == &(ds->out[part_idx]->partAdList->neighbors[sub_idx][size_out_neigh - 1]))
//				cursor = nullptr;
//			else
//				cursor = cursor + 1;
//		}
//
//		return *this;
//	}
//
//	NodeID operator*() {
//		return cursor->getNodeID();
//	}
//
//	Weight extractWeight() {
//		return cursor->getWeight();
//	}
//};



 // ------------------------------adList_chunk---------------------------------- NEW
template<typename U>
class neighborhood_iter<adListChunked<U>> {
	friend class neighborhood<adListChunked<U>> ;
private:
	U *cursor;
	U *endPtr;

public:

	neighborhood_iter(adListChunked<U> *ds, NodeID node, bool in_neigh) {
		cursor = nullptr;
		endPtr = nullptr;

		int64_t part_idx = ds->pt_hash(node);
		int64_t sub_idx = ds->hash_within_chunk(node);

		if(in_neigh){
			const int sz = ds->in[part_idx]->partAdList->neighbors[sub_idx].size();
			if(sz){
				cursor = &(ds->in[part_idx]->partAdList->neighbors[sub_idx][0]);
				endPtr = cursor + sz;
			}
		}
		else{
			const int sz =  ds->out[part_idx]->partAdList->neighbors[sub_idx].size();
			if(sz){
				cursor = &( ds->out[part_idx]->partAdList->neighbors[sub_idx][0]);
				endPtr = cursor + sz;
			}
		}

		while(cursor < endPtr){
			if(cursor->node >= 0){
				//found a valid node
				return;
			}
			cursor++;
		}
		cursor = nullptr;
		endPtr = nullptr;
	}

	bool operator!=(const neighborhood_iter<adListChunked<U>> &it) {
		return cursor != it.cursor;
	}

	neighborhood_iter& operator++() {
		cursor++;
		while(cursor < endPtr){
			if(cursor->node >= 0){
				//found a valid node
				return *this;
			}
			cursor++;
		}
		cursor = nullptr;
		return *this;
	}

	neighborhood_iter& operator++(int) {
		cursor++;
		while(cursor < endPtr){
			if(cursor->node >= 0){
				//found a valid node
				return *this;
			}
			cursor++;
		}
		cursor = nullptr;
		return *this;
	}

	NodeID operator*() {
		return cursor->getNodeID();
	}

	Weight extractWeight() {
		return cursor->getWeight();
	}
};


template<typename U>
class neighborhood_iter<darhh<U>> {
	friend class neighborhood<darhh<U>> ;
private:
	hd_rhh<U> *hd;
	ld_rhh<U> *ld;
	typename hd_rhh<U>::iter hd_iter;
	typename ld_rhh<U>::iter ld_iter;
	bool low_degree;
public:
	inline neighborhood_iter& operator=(neighborhood_iter const &it);
	inline bool operator!=(neighborhood_iter const &it);
	neighborhood_iter& operator++();
	neighborhood_iter& operator++(int);
	inline NodeID operator*();
	inline Weight extractWeight();
	void set_begin(darhh<U> *ds, NodeID n, bool in);
	void set_end();
};

template<typename U>
void neighborhood_iter<darhh<U>>::set_begin(darhh<U> *ds, NodeID src, bool in) {
	if (in) {
		ld = ds->in[ds->pt_hash(src)]->ld;
		hd = ds->in[ds->pt_hash(src)]->hd;
	} else {
		ld = ds->out[ds->pt_hash(src)]->ld;
		hd = ds->out[ds->pt_hash(src)]->hd;
	}
	low_degree = ld->get_degree(src);
	if (low_degree)
		ld_iter = ld->begin(src);
	else
		hd_iter = hd->begin(src);
}

template<typename U>
void neighborhood_iter<darhh<U>>::set_end() {
	ld_iter.cursor = nullptr;
	hd_iter.cursor = nullptr;
}

template<typename U>
neighborhood_iter<darhh<U>>&
neighborhood_iter<darhh<U>>::operator=(neighborhood_iter const &other) {
	hd = other.hd;
	ld = other.ld;
	hd_iter = other.hd_iter;
	ld_iter = other.ld_iter;
	low_degree = other.low_degree;
	return *this;
}

template<typename U>
bool neighborhood_iter<darhh<U>>::operator!=(neighborhood_iter const &it) {
	if (low_degree)
		return ld_iter != it.ld_iter;
	else
		return hd_iter != it.hd_iter;
}

template<typename U>
neighborhood_iter<darhh<U>>& neighborhood_iter<darhh<U>>::operator++() {
	if (low_degree)
		++ld_iter;
	else
		++hd_iter;
	return *this;
}

template<typename U>
neighborhood_iter<darhh<U>>& neighborhood_iter<darhh<U>>::operator++(int) {
	if (low_degree)
		++ld_iter;
	else
		++hd_iter;
	return *this;
}

template<typename U>
NodeID neighborhood_iter<darhh<U>>::operator*() {
	if (low_degree)
		return ld_iter.cursor->getNodeID();
	else
		return hd_iter.cursor->getNodeID();
}

template<typename U>
Weight neighborhood_iter<darhh<U>>::extractWeight() {
	if (low_degree)
		return ld_iter.cursor->getWeight();
	else
		return hd_iter.cursor->getWeight();
}

template<typename T>
class neighborhood {
private:
	T *ds;
	NodeID node;
	bool in_neigh;
public:
	neighborhood(NodeID _node, T *_ds, bool _in_neigh) :
			ds(_ds), node(_node), in_neigh(_in_neigh) {
	}
	neighborhood_iter<T> begin() {
		return neighborhood_iter<T>(ds, node, in_neigh);
	}
	neighborhood_iter<T> end() {
		neighborhood_iter<T> n = neighborhood_iter<T>(ds, node, in_neigh);
		n.cursor = nullptr;
		return n;
	}
};

//#if defined(USE_HYBRID_HASHMAP_WITH_GROUPING)								\
//	|| defined(USE_HYBRID_HASHMAP_WITH_GROUPING_AND_EDGE_ARR_LOCKING)		\
//	|| defined(USE_HYBRID_HASHMAP_WITH_GROUPING_TIGHTER)

#if defined(USE_CAHCE_FRIENDLY_HASH_ONLY)

template<typename U>
class neighborhood<GraphTango<U>> {
private:
	U* _start;
	uint64_t _size;
public:
	neighborhood(NodeID _node, GraphTango<U> *_ds, bool _in_neigh) {
		Vertex<U> &v = _ds->vArray[_node];
		if(_in_neigh){
			GraphTangoHash<U> &edges = v.inEdges;
			_size = edges.degree;
			if(edges.adjList == nullptr){
				//create adjacency list from hash table
				edges.adjList = (U*)globalAllocator.allocate(_size * sizeof(U));
				u32 idx = 0;
				const u32 cap = edges.capacity;
				for(u32 j = 0; j < cap; j++){
					if(edges.neighArr[j].node < FLAG_TOMB_STONE){
						edges.adjList[idx++] = edges.neighArr[j];
					}
				}
				//assert(idx == _size);
			}
			_start = edges.adjList;
		}
		else{
			GraphTangoHash<U> &edges = v.outEdges;
			_size = edges.degree;
			if(edges.adjList == nullptr){
				//create adjacency list from hash table
				edges.adjList = (U*)globalAllocator.allocate(_size * sizeof(U));
				u32 idx = 0;
				const u32 cap = edges.capacity;
				for(u32 j = 0; j < cap; j++){
					if(edges.neighArr[j].node < FLAG_TOMB_STONE){
						edges.adjList[idx++] = edges.neighArr[j];
					}
				}
				//assert(idx == _size);
			}
			_start = edges.adjList;
		}
	}
	neighborhood_iter<GraphTango<U>> begin() {
		return neighborhood_iter<GraphTango<U>>(_start);
	}
	neighborhood_iter<GraphTango<U>> end() {
		return neighborhood_iter<GraphTango<U>>(_start + _size);
	}
};

#elif defined(USE_HYBRID_HASHMAP_WITH_GROUPING)								\
	|| defined(USE_HYBRID_HASHMAP_WITH_GROUPING_AND_EDGE_ARR_LOCKING)		\
	|| defined(USE_HYBRID_HASHMAP_WITH_GROUPING_TIGHTER)					\
	|| defined(USE_GT_BALANCED_TYPE3_ONLY)

template<typename U>
class neighborhood<GraphTango<U>> {
private:
	U* _start;
	uint64_t _size;
public:
	neighborhood(NodeID _node, GraphTango<U> *_ds, bool _in_neigh) {
		if(_in_neigh){
			_start = _ds->vArray[_node].inEdges.neighArr;
			_size = _ds->vArray[_node].inEdges.degree;
		}
		else{
			_start = _ds->vArray[_node].outEdges.neighArr;
			_size = _ds->vArray[_node].outEdges.degree;
		}
	}
	neighborhood_iter<GraphTango<U>> begin() {
		return neighborhood_iter<GraphTango<U>>(_start);
	}
	neighborhood_iter<GraphTango<U>> end() {
		return neighborhood_iter<GraphTango<U>>(_start + _size);
	}
};

#elif 	defined(USE_GT_BALANCED)					\
		|| defined(USE_GT_BALANCED_MALLOC) 			\
		|| defined(USE_GT_BALANCED_STDMAP) 			\
		|| defined(USE_GT_BALANCED_MALLOC_STDMAP)	\
		|| defined(USE_GT_BALANCED_DYN_PARTITION)	\
		|| defined(USE_GT_BALANCED_ABSEIL)			\
		|| defined(USE_GT_BALANCED_RHH)				\
		|| defined(USE_GT_BALANCED_TSL_RHH)

template<typename U>
class neighborhood<GraphTango<U>> {
private:

	uint64_t degree;
	U* arr = nullptr;
	U** blocks = nullptr;
	bool isType3 = false;


public:
	neighborhood(NodeID _node, GraphTango<U> *_ds, bool _in_neigh) {
		if(_in_neigh){
			if(_ds->vArray[_node].inEdges.capacity <= EdgeArray<U>::TH0){
				arr = _ds->vArray[_node].inEdges.etype.type1.neigh;
			}
			else if(_ds->vArray[_node].inEdges.capacity <= EdgeArray<U>::TH1){
				arr = _ds->vArray[_node].inEdges.etype.type2.neighArr;
			}
			else{
				arr = _ds->vArray[_node].inEdges.etype.type3.blockList[0];
				blocks = _ds->vArray[_node].inEdges.etype.type3.blockList.data();
				isType3 = true;
			}
			degree = _ds->vArray[_node].inEdges.degree;
		}
		else{
			if(_ds->vArray[_node].outEdges.capacity <= EdgeArray<U>::TH0){
				arr = _ds->vArray[_node].outEdges.etype.type1.neigh;
			}
			else if(_ds->vArray[_node].outEdges.capacity <= EdgeArray<U>::TH1){
				arr = _ds->vArray[_node].outEdges.etype.type2.neighArr;
			}
			else{
				arr = _ds->vArray[_node].outEdges.etype.type3.blockList[0];
				blocks = _ds->vArray[_node].outEdges.etype.type3.blockList.data();
				isType3 = true;
			}
			degree = _ds->vArray[_node].outEdges.degree;
		}
	}
	neighborhood_iter<GraphTango<U>> begin() const {
		return neighborhood_iter<GraphTango<U>>(0, arr, blocks, isType3);
	}
	neighborhood_iter<GraphTango<U>> end() const {
		return neighborhood_iter<GraphTango<U>>(degree, arr, blocks, isType3);
	}
};

#elif defined(USE_GT_UPDATE)


#else

template<typename U>
class neighborhood<GraphTango<U>> {
private:
	U* _start;
	uint64_t _size;
public:
	neighborhood(NodeID _node, GraphTango<U> *_ds, bool _in_neigh) {
		if(_in_neigh){
			_start = _ds->vArray.inNeighArr[_node];
			_size = _ds->vArray.inDegree[_node];
		}
		else{
			_start = _ds->vArray.outNeighArr[_node];
			_size = _ds->vArray.outDegree[_node];
		}
	}
	neighborhood_iter<GraphTango<U>> begin() {
//		#pragma omp atomic
//		g_totalEdges += _size;
		return neighborhood_iter<GraphTango<U>>(_start);
	}
	neighborhood_iter<GraphTango<U>> end() {
		return neighborhood_iter<GraphTango<U>>(_start + _size);
	}
};

#endif


template<typename U>
class neighborhood<darhh<U>> {
private:
	using iter = neighborhood_iter<darhh<U>>;
	NodeID src;
	darhh<U> *ds;
	bool in;
public:
	neighborhood(NodeID src, darhh<U> *ds, bool in) :
			src(src), ds(ds), in(in) {
	}
	iter begin() {
		iter it;
		it.set_begin(ds, src, in);
		return it;
	}
	iter end() {
		iter it;
		it.set_end();
		return it;
	}
};

template<typename T>
neighborhood<T> in_neigh(NodeID n, T *ds) {
	if (ds->directed)
		return neighborhood<T>(n, ds, true);
	else
		return neighborhood<T>(n, ds, false);
}

template<typename T>
neighborhood<T> out_neigh(NodeID n, T *ds) {
	return neighborhood<T>(n, ds, false);
}

#endif // TRAVERSAL_H_

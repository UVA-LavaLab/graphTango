#pragma once

#include "common.h"
#include "CustomAllocator.h"

template <typename Key, typename Value>
class CacheFriendlyHash{

	typedef struct {
		Key key;
		Value val;
	} KeyValuePair;

	u32 degree = 0;
	u32 capacity = 0;
	KeyValuePair* __restrict data = nullptr;



	void insert(const Key& key, const Value& value){
		//check capacity

	}

	Value& getValue(const Key& key){

	}

	const Value& getValue(const Key& key) const {

	}



};

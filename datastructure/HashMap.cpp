/*
 * This is a data structure similar to hashmap . This map can hold the size number of matrix having key k. The main purpose of the class is to serve as a caching.
 * @author Paras Babu Tiwari
 * HashMap.cpp
 *
 *  Created on: Jan 15, 2012
 *      Author: user
 */

#include "HashMap.h"

/**
 * Create a hashmap.
 * @param size: size of the hashmap
 */
HashMap::HashMap(int size) {
	// TODO Auto-generated constructor stub
	id = new int[size];
	marray = new Matrices*[size];
	this->s = size;
	current = -1;

}

HashMap::~HashMap() {
	// TODO Auto-generated destructor stub
	delete id;
	delete marray;
}
/*
 * Add an entry in the hashmap.
 * We maintained a circular list of id, and the map contains last size items.
 */
void HashMap::add(int key, Matrices *value) {
	if (key > 0) {
		current++;
		current = current % s;
		assert(current<s);
		id[current] = key;
		marray[current] = value;
	}

}
/*
 * Get an item.
 */
Matrices * HashMap::get(int key) {
	if(key<0)
		return NULL;
	for (int i = 0; i <= current; i++) {
		if (id[i] == key)
			return marray[i];
	}
	return NULL;
}

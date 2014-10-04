/*
 * HashMap.h
 *
 *  Created on: Jan 15, 2012
 *      Author: user
 */

#ifndef HASHMAP_H_
#define HASHMAP_H_
#include "Matrices.h"
class HashMap {
public:
	HashMap(int );
	virtual ~HashMap();
	int *id;
	Matrices **marray;
	int current;
	int s;

	void add(int key,Matrices *value);
	Matrices* get(int key);


};

#endif /* HASHMAP_H_ */

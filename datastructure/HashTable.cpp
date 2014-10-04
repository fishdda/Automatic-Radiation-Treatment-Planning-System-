/** @file HashTable.cpp
 *
 *  @author Jon Turner
 *  @date 2011
 *  This is open source software licensed under the Apache 2.0 license.
 *  See http://www.apache.org/licenses/LICENSE-2.0 for details.
 */

#include "HashTable.h"

/** Constructor for HashTable, allocates space and initializes table.
 *  N1 is the limit on the range of values; it must be less than 2^20.
 */
HashTable::HashTable() {
	n = MAXSIZE;

	for (nb = 1; 8 * nb <= n; nb <<= 1) {
	}
	nb = max(nb, 4);
	bktMsk = nb - 1;
	kvxMsk = (8 * nb) - 1;
	fpMsk = ~kvxMsk;

	bkt = new bkt_t[2 * nb];
	pairs = new KeyValPairs[n + 1];
	kvx = new UiSetPair(n);

	clear();
}
;

/** Destructor for HashTable. */
HashTable::~HashTable() {
	delete[] bkt;
	delete[] pairs;
	delete kvx;
}

/** Clear the hashtable contents. */
void HashTable::clear() {
	for (int i = 0; i < 2 * nb; i++) {
		for (int j = 0; j < BKT_SIZ; j++)
			bkt[i][j] = 0;
	}
	kvx->reset();
}

/** Compute a bucket index and fingerprint, for a given key.
 *
 *  Hashit uses multiplicative hashing with one of two
 *  different multipliers, after first converting
 *  the 64 bit integer into a 32 bit integer.
 *
 *  @param key is the key to be hashed
 *  @param hf is either 0 or 1 and selects one of two hash functions
 *  @param b is a reference; on return its value is equal to the
 *  hash bucket for the given key
 *  @param fp is a reference; on return its value is equal to the
 *  fingerprint for the given key
 */
void HashTable::hashit(uint64_t key, int hf, uint32_t& b, uint32_t& fp) {
	const uint32_t A0 = 0xa96347c5;
	const uint32_t A1 = 0xe65ac2d3;

	uint32_t x, y;
	uint64_t z;

	x = (key >> 16) & 0xffff0000;
	x |= (key & 0xffff);
	y = (key >> 48) & 0xffff;
	y |= (key & 0xffff0000);
	z = x ^ y;
	z *= (hf == 0 ? A0 : A1);
	b = (z >> 16) & bktMsk;
	fp = (z >> 13) & fpMsk;
}

/** Get the value for a specified key.
 *  @param key is the keey to be looked up in the table
 *  @return the value stored for the given key, or UNDEF_VAL
 *  if there is none.
 */
double HashTable::get(uint64_t key) {
	uint32_t b, kvIndex, fp;

	// check bucket in the first half of the bucket array
	hashit(key, 0, b, fp);

	for (int i = 0; i < BKT_SIZ; i++) {
		if ((bkt[b][i] & fpMsk) == fp) {
			kvIndex = bkt[b][i] & kvxMsk;

			if (pairs[kvIndex].key == key) {
				return pairs[kvIndex].val;
			}
		}
	}

	// check bucket in the second half of the bucket array
	hashit(key, 1, b, fp);

	b += nb;

	for (int i = 0; i < BKT_SIZ; i++) {
		if ((bkt[b][i] & fpMsk) == fp) {
			kvIndex = bkt[b][i] & kvxMsk;
			if (pairs[kvIndex].key == key) {
				return pairs[kvIndex].val;
			}
		}
	}
	return UNDEF_VAL;
}

/** Put a (key,value) pair into the map.
 *  If the key matches one that is already in the map,
 *  the previously associated value is replaced.
 *  @param key is the key part of the pair
 *  @param val is the value part of the pair
 *  @return true on success, false on failure.
 */
bool HashTable::put(uint64_t key, double val) {
	int i, j0, j1, n0, n1;
	uint32_t b0, b1, fp0, fp1;

	// Count the number of unused items in each bucket
	// and find an unused item in each (if there is one)
	hashit(key, 0, b0, fp0);
	n0 = 0;
	for (i = 0; i < BKT_SIZ; i++) {
		if (bkt[b0][i] == 0) {
			n0++;
			j0 = i;
		} else if ((bkt[b0][i] & fpMsk) == fp0) {
			// key may already in map
			int kvIndex = bkt[b0][i] & kvxMsk;
			if (pairs[kvIndex].key == key) {
				// key is in map, update value
				pairs[kvIndex].val = val;
				return true;
			}
		}
	}
	hashit(key, 1, b1, fp1);
	b1 += nb;
	n1 = 0;
	for (i = 0; i < BKT_SIZ; i++) {
		if (bkt[b1][i] == 0) {
			n1++;
			j1 = i;
		} else if ((bkt[b1][i] & fpMsk) == fp1) {
			// key may already in map
			int kvIndex = bkt[b1][i] & kvxMsk;
			if (pairs[kvIndex].key == key) {
				// key is in map, update value
				pairs[kvIndex].val = val;
				return true;
			}
		}
	}
	// If no unused entry in either bucket, give up.
	if (n0 + n1 == 0) {
		printf("No un-used entry\n");
		return false;
	}

	// store the key value in keyVec and add entry in least-loaded bucket
	int kvIndex = kvx->firstOut();
	if (kvIndex == 0) {
		printf("No more room\n");
		exit(0);
		return false; // no more room
	}
	kvx->swap(kvIndex);
	pairs[kvIndex].key = key;
	pairs[kvIndex].val = val;
	if (n0 >= n1)
		bkt[b0][j0] = fp0 | (kvIndex & kvxMsk);
	else
		bkt[b1][j1] = fp1 | (kvIndex & kvxMsk);

	return true;
}
/*
 * Initialize the Data structure.
 * @Param
 * 	r: List of row.
 * 	c: List of col.
 * 	val: List of value of the matrix. The sparse matrix has a non-zero entry at r(k),c(k) and its value is val(k).
 * 	nrow: # of non-zero element.
 * 	@Author: Paras Babu Tiwari
 */
void HashTable::initialize(double *r, double *c, double *val, int nrow) {
	for (int i = 0; i < nrow; i++) {
		insertItem(r[i], c[i] - 1, val[i]);
	}
}
/*
 * Insert the item in the hash.
 * @Param:
 * 	i: Row #.
 * 	j: Column #.
 * 	val: Value of (i,j) entry.
 *
 */
void HashTable::insertItem(int i, int j, double val) {

	uint64_t key = makeKey(i, j);
	this->put(key, val);
}
uint64_t HashTable::makeKey(uint64_t i, int j) {

	uint64_t key;
	string c; // string which will contain the result

	ostringstream cc; // stream used for the conversion
	cc << j;
	c = cc.str();
	int s = c.size();
	s = s + 1;
	int m = pow(10.0, s);
	//printf("r:%llu\n",r);
	i = i * m;
	//printf("m:%d\n",m);
	//printf("r:%llu\n",r);
	key = i + j;
	return key;
}
double HashTable::search(int r, int c) {

	uint64_t key = makeKey(r, c);

	double val = this->get(key);
	if (val == UNDEF_VAL)
		val = 0;
//	if(val != 0)
//		printf("row:%d col:%d val:%lf\n",r,c,val);
	return val;

}

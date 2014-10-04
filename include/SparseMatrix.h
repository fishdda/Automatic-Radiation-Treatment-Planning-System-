/*
 * SparseMatrix.h
 *
 *  Created on: Nov 20, 2011
 *      Author: user
 */

#include <stdlib.h>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include "Matrices.h"

#ifndef SPARSEMATRIX_H_
#define SPARSEMATRIX_H_
namespace ublas = boost::numeric::ublas;

class SparseMatrix {
	ublas::compressed_matrix<double, ublas::row_major> *value;
	//ublas::compressed_matrix<double, ublas::column_major> *value;
	 std::map<int,Matrices*> mp;
	//Matrices **indvMatrix;
	//double *value;
	//int row;
	//int col;

	long unsigned int nzEntry;	//Holds the number of non-zero entry
	//int numRow,numCol;
	//int numFullRow,numFullCol;
//	HashTable *htbl;
//	HashMap *hmap;
	double *rowIndex,*colIndex;
public:

	SparseMatrix(double* r, double* c, double *val,long unsigned int nonzeroel,
			long unsigned int rf, long unsigned int cf);
	//SparseMatrix(double* row,double* col,double *value,int nrow,int rf,int cf);
	virtual ~SparseMatrix();
	Matrices* getSubMatrix(int id, double *rows,long unsigned int numRow);

	long unsigned int getRow() {
		return (*value).size1();
	}
	long unsigned int getCol() {
		return (*value).size2();
	}

	//void printMatrix();
	double getAt(long unsigned int i, long unsigned int j);
	//int searchr(int r,int c);
	long unsigned int getNzEntry() {
		return nzEntry;
	}
	double getValAt(long unsigned int i);
	unsigned long int getRowNumber(long unsigned int i);
	unsigned long int getColNumber(long unsigned int i);

};

#endif /* SPARSEMATRIX_H_ */

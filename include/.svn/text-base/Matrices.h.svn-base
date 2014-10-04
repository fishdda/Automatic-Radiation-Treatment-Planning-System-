/*
 * DoubleMatrices.h
 *
 *  Created on: Nov 19, 2011
 *      Author: user
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <assert.h>
#include <vector>
#include "Vector.h"
#include "stdinc.h"
#ifndef MATRIX_H_
#define MATRIX_H_
#undef NDEBUG
namespace ublas = boost::numeric::ublas;
class Matrices {
	//double **mat;
	ublas::compressed_matrix<double, ublas::row_major>* mat;
	
	std::vector<long int> rowNumber;
	std::vector<long int> colNumber;
	
public:
	int row, col;
	std::vector<double> value;
	std::vector<double> negValue;
	Matrices();
	Matrices(ublas::compressed_matrix<double, ublas::row_major> *m);
	Matrices(const Matrices* cpm);
	//Matrices(double* r,double* c,double *val,int nrow,int ncol,int);
	//Matrices(boost::numeric::ublas::matrix_range < double > m);

	//Matrices(boost::numeric::ublas::matrix_range < boost::numeric::ublas::compressed_matrix<double> > m);
	//Matrices(const double *x,int r);
	//Matrices(const double *x,int r,int c);
	~Matrices();
	//double **getMat(){return mat;}
	void createMatrices();
	void print();
	void printNonZeroEntry();
	void copyNzEl(int *iRow, int *iCol, double *value,int &l, int &total);
	int getRow() {
		return row;
	}
	int getCol() {
		return col;
	}
	PrioritizedMain::Vector* matrixMultiply(PrioritizedMain::Vector *v);
	void copy(Matrices *other, int start);
	void transpose();
	double sum(int, int);
	void colSums(double *result, int &total,int n);
	void colSums(double *result);
	double sum();
	void divide(double factor);
	void multiply(double factor);
	PrioritizedMain::Vector* matrixTransMultiply(PrioritizedMain::Vector *v);
	void substract(double factor);
	void substract(double factor, int, int);
	void subfrom(double factor, int start, int end);
	void subfrom(double factor);
	void getDouble(double *x);
	void getSubMatrices(double *rows, int numRow, Matrices *result);
	double multiplyAt(int i, int j, double val);
	double getAt(int r, int c) {
		//printf("r:%d\n",r);
		//printf("c:%d\n",c);
		return (*mat)(r, c);
	}
	;
	double average(int s, int e);
	double getMax(int s, int e);
	double getMin(int s, int e);
	double average();
	double getMax();
	double getMin();
	void setAt(int row, int col, double val);
	double sum(boost::numeric::ublas::vector<double> v);
//	double average();
//		double getMax();
//		double getMin();
	void vectorMultiply(Matrices *mat1, Matrices *mat2, Matrices *result);
	void divide(double factor, int start, int end);
	PrioritizedMain::Vector* getVector(boost::numeric::ublas::vector<double> v);
	void addRowNumber(long int r) {
		rowNumber.push_back(r);
	}
	void addColNumber(long int c) {
		colNumber.push_back(c);
	}
	void addElement(double v)
	{
		value.push_back(v);
		
	}
	std::vector<long int> getNzRows()
	{
		return rowNumber;
	}
	std::vector<long int> getNzCols()
	{
		return colNumber;
	}
	std::vector<double> getNzElements()
	{
		return value;
	}
	long int getRowAt(int i) {
		return rowNumber[i];
	}
	long int getColAt(int i) {
		return colNumber[i];
	}
	long int getTotalNonZero() {
		return rowNumber.size();
	}
	void setMatrix(ublas::compressed_matrix<double, ublas::row_major> *m) {
		mat = m;
		row = m->size1();
		col = m->size2();
	}
	void printMatAdd();

};

#endif /* MATRIX_H_ */

/*
 * Vector.h
 *
 *  Created on: Nov 28, 2012
 *      Author: parastiwari
 */

#ifndef VECTOR_H_
#define VECTOR_H_

#include <boost/numeric/ublas/vector.hpp>
#include <algorithm>
#include <stdio.h>
namespace PrioritizedMain {

class Vector {
public:
	Vector();
	Vector(int s);
	Vector(const double *x,const size_t s);
	Vector(boost::numeric::ublas::vector<double> ublasVector);
	virtual ~Vector();
	boost::numeric::ublas::vector<double> *v;
	double getMin(int s,int e);
	double getMin();
	double getMax(int s,int e);
	double getMax();
	long double sum();
	double average();
	double *getDblArray();
	inline double getAt(int i){return V()(i);}
	int size();
	Vector* operator- (double factor);
	Vector* operator/(double factor);
	double operator*(Vector *vec);
	Vector *square();
	Vector* operator*(double factor);
	Vector* subFrom(double factor);
	inline boost::numeric::ublas::vector<double> V(){return *v;}
	void print();
	void printAll();
	void sort();
};

} /* namespace PrioritizedMain */
#endif /* VECTOR_H_ */

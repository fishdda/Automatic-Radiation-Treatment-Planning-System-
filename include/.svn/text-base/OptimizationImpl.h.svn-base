/*
 * OptFormImpl.h
 *
 *  Created on: Jul 11, 2012
 *      Author: user
 */
#include "ProblemReps.h"
#include "Matrices.h"
#ifndef OPTFORMIMPL_H_
#define OPTFORMIMPL_H_
#undef NDEBUG
#include <assert.h>

class OptimizationImpl  {
	int step;
	ProblemReps *current;
public:
	OptimizationImpl(int step,ProblemReps *prob);
	virtual ~OptimizationImpl();
	long double evaluateObjective(const double *x);
	void evaluategradf(const double *x, double *g);
	void evalg(const double *x,double *g);
	int getNumVar();
	int getNumConstraint();
	int getNumJac();
	bool getHessian(int *iRow, int *iCol, double *value, double obj_factor,const double* lambda);
	int getNumHess();
	void getJacobian(int *, int *, double*,const double *);
	void getStartingPoint(double *x,double *z_l,double *z_u);
	void getBounds(double *x_l, double *x_u,double *g_l,double *g_u);
	void PrintDose();
	void setW(const double *x,const double *z_l,const double *z_u,const double *lambda);
	void setObjectiveFunValue(double val);
	//void setStep(int s,ProblemReps *prob);
	//void printDose(Matrices *w);
	//int setw();

};

#endif /* OPTFORMIMPL_H_ */

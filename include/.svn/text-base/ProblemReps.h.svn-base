/*
 * Abstract class representing the problem.
 * This class can't be instantiate and must be override and implement the functions representing the problem.
 * ProblemReps.h
 *
 *  Created on: Jul 16, 2012
 *      Author: Paras Babu Tiwari
 */

#ifndef PROBLEMREPS_H_
#define PROBLEMREPS_H_

class ProblemReps {

public:
	ProblemReps();
	virtual ~ProblemReps();
	virtual long double evaluateObjective(const double *x) = 0;
	virtual void evaluategradf(const double *x, double *g)=0;
	virtual bool getHessian(int *iRow, int *iCol, double *value,
			double obj_factor,const double *lambda);
	virtual void getJacobian(int *, int *, double*,const double* )=0;

	virtual void evalg(const double *x, double *g)=0;
	virtual int getNumVar() =0;
	virtual int getNumConstraint()=0;
	virtual int getNumJac() = 0;
	virtual int getNumHess();
	virtual void getBounds(double *x_l, double *x_u, double *g_l,
			double *g_u)=0;
	virtual void printDose()=0;
	virtual void setW(const double *x,const double *z_l,const double *z_u,const double *lambda)=0;
	virtual void setObjectiveFunValue(double val)=0;
	virtual void getStartingPoint(double *x,double *z_l,double *z_u);
};

#endif /* PROBLEMREPS_H_ */

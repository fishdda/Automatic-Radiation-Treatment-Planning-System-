/*
 * StepIIImpl.h
 *
 *  Created on: Sep 15, 2012
 *      Author: user
 */

#ifndef STEPIIIMPL_H_
#define STEPIIIMPL_H_
#include "ProblemReps.h"
#include "StepIImpl.h"
#include "Matrices.h"
#include "stdinc.h"
#ifdef dummy
#include "oneVoxel.h"
#else
#include "MatlabCommunicationManager.h"
#endif

class StepIIImpl: public ProblemReps {
#ifdef dummy
	oneVoxel *mtlb;
#else
	MatlabCommunicationManager *mtlb;
#endif
	//oneVoxel *mtlb;
	double *grad;
	int NUMBEAMLET;
	StepIImpl *step2;
	double obj;
	PrioritizedMain::Vector *w3;
public:
	double *xn;
#ifdef dummy
	StepIIImpl(oneVoxel *mtlb, StepIImpl *stpp2);
#else
	StepIIImpl(MatlabCommunicationManager *m, StepIImpl *stpp2);
#endif
	//
	virtual ~StepIIImpl();
	virtual void getStartingPoint(double *x,double *z_l,double *z_u);
	virtual long double evaluateObjective(const double *x);
	virtual void evaluategradf(const double *x, double *g);
	virtual void getJacobian(int *, int *, double*, const double*);
	virtual void evalg(const double *x, double *g);
	void evalg(const double *x, double *g, int& total);
	void evalNewConst(const double *x, double *g, int& total);
	virtual int getNumVar();
	virtual int getNumConstraint();
	virtual int getNumJac();
	virtual int getNumHess();
	virtual void getBounds(double *x_l, double *x_u, double *g_l, double *g_u);
	void getBounds(double *x_l, double *x_u, double *g_l, double *g_u,
			double slip,double s2=0);
	void getJacobian(int *iRow, int *iCol, double *value, const double *x,
			int& total, int& l);
	virtual void printDose();
	void printDose(PrioritizedMain::Vector *w);
	virtual void setW(const double *x,const double *z_l,const double *z_u,const double *lambda);
	virtual void setObjectiveFunValue(double val) {
		obj = val;
	}
	bool getHessian(int *iRow, int *iCol, double *value, double obj_factor,
			const double *lambda);
	bool getHessian(int *iRow, int *iCol, double *value,
			double obj_factor, const double *lambda,int& total,bool step4);
	void handleNewConst(int *iRow, int *iCol, double *value, const double *x,
			int &l, int &total);
	void divideByVoxel(PrioritizedMain::Vector **D);
	int getNumNewConst();
	int getNumNewConstJac();
	void divideByVoxel(Matrices **m);
	PrioritizedMain::Vector *getW3(){return w3;}
	StepIImpl* Step2(){return step2;}
};

#endif /* STEPIIIMPL_H_ */

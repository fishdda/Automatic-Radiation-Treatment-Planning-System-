/*
 * StepIVmpl.h
 *
 *  Created on: Sep 15, 2012
 *      Author: user
 */

#ifndef StepIVmpl_H_
#define StepIVmpl_H_
#include "ProblemReps.h"
#include "StepIIImpl.h"
#include "Matrices.h"
#include "stdinc.h"
#ifdef dummy
	#include "oneVoxel.h"
#else
	#include "MatlabCommunicationManager.h"
#endif

class StepIVmpl: public ProblemReps {
#ifdef dummy
	oneVoxel *mtlb;
#else
	MatlabCommunicationManager *mtlb;
#endif
	//oneVoxel *mtlb;
	int NUMBEAMLET;
	StepIIImpl *step3;
	double obj;
	PrioritizedMain::Vector *w4;
	int numJacFromNewConst;
	double *newconsjac;
public:
#ifdef dummy
	StepIVmpl(oneVoxel *mtlb, StepIIImpl *stpp2);
#else
	StepIVmpl(MatlabCommunicationManager *m, StepIIImpl *stpp2);
#endif
	//
	virtual ~StepIVmpl();
	virtual void getStartingPoint(double *x,double *z_l,double *z_u);
	virtual long double evaluateObjective(const double *x);
	virtual void evaluategradf(const double *x, double *g);
	virtual void getJacobian(int *, int *, double*, const double*);
	virtual void evalg(const double *x, double *g);
	void evalNewConst(const double *x, double *g,int& total);
	virtual int getNumVar();
	virtual int getNumConstraint();
	virtual int getNumJac();
	virtual int getNumHess();
	virtual void getBounds(double *x_l, double *x_u, double *g_l, double *g_u);
	virtual void printDose();
	virtual void setW(const double *x,const double *z_l,const double *z_u,const double *lambda);
	virtual void  setObjectiveFunValue(double val){obj = val;}
	bool getHessian(int *iRow, int *iCol, double *value, double obj_factor,
			const double *lambda);
	void handleNewConst(int *iRow, int *iCol, double *value,
			const double *x, int &l, int &total);

	int getNumNewConst();
	int getNumNewConstJac();
	void divideByVoxel(Matrices **m);

};

#endif /* StepIVmpl_H_ */

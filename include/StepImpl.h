/*
 * StepImpl.h
 *
 *  Created on: Nov 20, 2011
 *      Author: user
 */
#ifndef StepImpl_H_
#define StepImpl_H_

#include <algorithm>
#include "Matrices.h"
#include "stdinc.h"

#ifdef dummy
	#include "oneVoxel.h"
#else
	#include "MatlabCommunicationManager.h"
#endif
#include "ProblemReps.h"
//#include "oneVoxel.h"
#include "FileManager.h"
#include "ProblemReps.h"
#include <assert.h>
#include <stdio.h>
#include <iostream>

#undef NDEBUG
#include <assert.h>

class StepImpl: public ProblemReps {
#ifdef dummy
	//SmartPtr<oneVoxel>  mtlb;
	
	oneVoxel  *mtlb;
#else
	MatlabCommunicationManager *mtlb;
#endif
	//oneVoxel *mtlb;
	double *GS;
	PrioritizedMain::Vector *w1;
	int NUMBEAMLET;
	double obj;
	double *minDose,*maxDose;
	time_t fstart_t, fend_t,gstart_t,gend_t,jacstart_t,jacend_t,hessstart_t,hessend_t;
	Matrices **targetInf;
	int numTargets;
	double *ZL;
	double *ZU;
	int numPBs;
	std::vector<double> negValue;
public:
#ifdef dummy
	StepImpl(oneVoxel *m);
#else
	StepImpl(MatlabCommunicationManager *mtlb);
#endif
	virtual ~StepImpl();

	void copy(double *D, double **src, int t);
	void getTargetDose(const double *x, PrioritizedMain::Vector **D);
	void getTargetDose(PrioritizedMain::Vector *w, PrioritizedMain::Vector **D);
	//int getRow(int i);
	//Matrices *getwMatrix(const double *x);
	//int getRowStepI(int i);

	void print();
	void getStepIDose(Matrices *w,Matrices *D);
	virtual void getStartingPoint(double *x,double *z_l,double *z_u);
	//int getNumVoxels(int i);

	void divideByVoxel(PrioritizedMain::Vector **D);

	void setStep(int s);

	//void sumByTarget(Matrices *D, double *G);
	//void evaluateG(double*& G,PrioritizedMain::Vector *w);
	void evaluateG(double* G,PrioritizedMain::Vector *w);

	//void setw(Matrices *w){this->w = w;}
	double getMinDose(int tar);
	double getMaxDose(int tar);

	virtual long double evaluateObjective(const double *x);
	virtual void evaluategradf(const double *x, double *g);
	void evaluategradf1(const double *x, double *g);
	virtual bool getHessian(int *iRow, int *iCol, double *value,
			double obj_factor,const double *lambda);
	virtual void getJacobian(int *, int *, double *, const double*);

	virtual void evalg(const double *x, double *g);
	virtual int getNumVar();
	virtual int getNumConstraint();
	virtual int getNumJac();
	virtual int getNumHess();
	virtual void getBounds(double *x_l, double *x_u, double *g_l, double *g_u);
	virtual void printDose();
	void printDose(PrioritizedMain::Vector *w);
	void setW(const double *x,const double *z_l,const double *z_u,const double *lambda);
	virtual void setObjectiveFunValue(double val){obj = val;}
	PrioritizedMain::Vector *getW1(){return w1;}
	inline double *getZL(){return ZL;}
	inline double *getZU(){return ZU;}
	void subPrescDose(PrioritizedMain::Vector **D);
};

#endif /* OPTFORMIMPL_H_ */

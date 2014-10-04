/*
 * StepIImpl.h
 *
 *  Created on: Nov 20, 2011
 *      Author: user
 */
#ifndef StepIImpl_H_
#define StepIImpl_H_

#include <stdio.h>
#include <iostream>
#include "Matrices.h"
#include "ProblemReps.h"
#include "StepImpl.h"
#ifdef dummy
	#include "oneVoxel.h"
#else
	#include "MatlabCommunicationManager.h"
#endif
//#include "DummyMatlabComm.h"

#include "FileManager.h"
#include "Constant.h"
#undef NDEBUG
#include <assert.h>

class StepIImpl: public ProblemReps {

#ifdef dummy
	oneVoxel *mtlb;
#else
	MatlabCommunicationManager *mtlb;
#endif
	//oneVoxel *mtlb;
	StepImpl *step1;
	int ZINDEX;
	int PINDEX;
	int YINDEX;
	int NUMYS;
	int NUMPS;
	int NUMZS;
	int NUMBEAMLET;
	//int cache; //This variable is use to decide whether we should calculate the new Jacobian or use old one.
	int *irow;
	int *icol;
	double *ivalue;
	PrioritizedMain::Vector *w2;
	double obj;
	int numFifth;
	double miAlpha;
	double **fourthconstjac;
	int *fjacsz;
	int *fourthnumrow;
	int numTargets;
	Matrices **targetInf;
	
public:
	double *xn;
#ifdef dummy
	StepIImpl(oneVoxel *mtlb, StepImpl *stp1);
#else
	StepIImpl(MatlabCommunicationManager *m, StepImpl *stp1);
#endif
	int yindex(){return YINDEX;}
	int zindex(){return ZINDEX;}
	int pindex(){return PINDEX;}
	//StepIImpl(oneVoxel *m, StepImpl *stp1);
	virtual ~StepIImpl();

	void handleFirstConst(int *iRow, int *iCol, double *value, const double *x,
			int &, int &);
	void handleSecondConst(int *iRow, int *iCol, double *value, const double *x,
			int &, int &);
	void handleThirdConst(int *iRow, int *iCol, double *value, const double *x,
			int &, int &);
	void handleFourthConst(int *iRow, int *iCol, double *value, const double *x,
			int &, int &);
	void handleFifthConst(int *iRow, int *iCol, double *value, const double *x,
			int &, int &);

	void setW(const double *x,const double *z_l,const double *z_u,const double *lambda); 
	virtual void getStartingPoint(double *x,double *z_l,double *z_u);
	virtual long double evaluateObjective(const double *x);
	virtual void evaluategradf(const double *x, double *g);

	virtual void getJacobian(int *, int *, double*, const double*);
	virtual bool getHessian(int *iRow, int *iCol, double *value,
			double obj_factor, const double *lambda);
	bool getHessian(int *iRow, int *iCol, double *value,
			double obj_factor, const double *lambda,int &total,bool step4);
	double getMialphaMax();

	virtual void evalg(const double *x, double *g);
	virtual int getNumVar();
	virtual int getNumConstraint();
	virtual int getNumHess();
	virtual void setObjectiveFunValue(double val){obj = val;}
	int getNumFirstConst();
	int getNumSecondConst();
	int getNumThirdConst();
	int getNumFourthConst();
	int getNumFifthConst();
	int getNumJacFromFirstConst();
	int getNumJacFromSecondConst();
	int getNumJacFromThirdConst();
	int getNumJacFromFourthConst();
	int getNumJacFromFifthConst();
	void eval1stConst(const double *x, int &total, double *g);
	void eval2ndConst(const double *x, int &total, double *g);
	void eval3rdConst(const double *x, int &total, double *g);
	void eval4thConst(const double *x, int &total, double *g);
	void eval5thConst(const double *x, int &total, double *g);
	void getJacobian(int *iRow, int *iCol, double *value, const double *x,
			int& total,int& l);
	void getBounds(double *x_l, double *x_u, double *g_l, double *g_u,double slip,double s2=0);
	void evalg(const double *x, double *g,int& total);
	virtual int getNumJac();
	virtual void printDose();

	virtual void getBounds(double *x_l, double *x_u, double *g_l, double *g_u);
	double getObjectiveFnValue(){return obj;}
	void printDose(PrioritizedMain::Vector *);
	StepImpl* Step1(){return step1;}
	PrioritizedMain::Vector *getW2(){return w2;}
};

#endif /* OPTFORMIMPL_H_ */

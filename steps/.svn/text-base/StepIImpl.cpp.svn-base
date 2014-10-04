/*
 * StepIImpl.cpp *
 *  Created on: Nov 20, 2011
 *      Author: Paras Babu Tiwari
 */

#include "StepIImpl.h"

#define ARTIFICIALVAR 3
#define FILESIZE 4
#ifdef dummy
StepIImpl::StepIImpl(oneVoxel *m, StepImpl *stp1)
#else
		StepIImpl::StepIImpl(MatlabCommunicationManager *m, StepImpl *stp1)
#endif
		{

	mtlb = m;
	numTargets = mtlb->getNumTargets();
	targetInf = new Matrices*[numTargets];
	printf("Space created\n");
	fjacsz = new int[numTargets];
	int *numvoxel = mtlb->getNumVoxel();
	for (int i = 0; i < numTargets; i++) {
		int id = mtlb->getTargetsId(i);
		int n = numvoxel[id];
		fjacsz[i] = 0;
		targetInf[i] = mtlb->influence->getSubMatrix(id, mtlb->voxels[id],
				numvoxel[mtlb->getTargetsId(i)]);
	}
	fourthnumrow = new int[numTargets];
	fourthconstjac = new double*[numTargets];
	this->step1 = stp1;
	//getColInfMat returns the # of beamlets.
	NUMBEAMLET = mtlb->getColInfMat();
	//We need to take into account of y's,p's and z's. I assumed that each organ has not more than one alpha value.So # y's equals to # of step2 organs.
	NUMYS = mtlb->getNumStepIIOar();
	//p's and z's equals to total number of voxels in step II organ * # of step II organs. This is because p and z has index over all voxel and step II organs.
	NUMPS = mtlb->getTotalStepIIVoxels() * mtlb->getNumStepIIOar();
	NUMZS = mtlb->getTotalStepIIVoxels() * mtlb->getNumStepIIOar();
	printf("NUMPS:%d\n", NUMPS);
	printf("NUMZS:%d\n", NUMZS);
	printf("NUMYS:%d\n", NUMYS);

	//+ mtlb->getTotalStepIIVoxels();
	//After beamlet index, I assumed y index
	YINDEX = NUMBEAMLET;
	//After yindex, I assumed pindex
	PINDEX = YINDEX + NUMYS;

	//After pindex, I assumed zindex
	ZINDEX = PINDEX + NUMPS; //mtlb->getColInfMat() + mtlb->getNumStepIIOar()
	printf("ZINDEX:%d\n", ZINDEX);
	printf("YINDEX:%d\n", YINDEX);
	printf("PINDEX:%d\n", PINDEX);

	irow = NULL;
	icol = NULL;
	ivalue = NULL;
	numFifth = 0;
	miAlpha = 0;
//	row = new int[getNumJac()];
//	col = new int[getNumJac()];
//	ivalue = new double[getNumJac()];
}
void StepIImpl::setW(const double *x,const double *z_l,const double *z_u,const double *lambda) {
		w2 = new PrioritizedMain::Vector(x, NUMBEAMLET);
		int nv = getNumVar();
		xn = new double[nv];
		memcpy(xn,x,sizeof(double)*nv);
		//w2 = new Matrices(x, NUMBEAMLET);
	}

void StepIImpl::getStartingPoint(double *x,double *z_l,double *z_u) {
	for (int i = 0; i < getNumVar(); i++) {
		x[i] = 0;
		
	}
	/*
	PrioritizedMain::Vector *w = step1->getW1();
	double *ZL = step1->getZL();
	double *ZU = step1->getZU();
	for (int i = 0; i < NUMBEAMLET; i++) {
		x[i] = w->getAt(i);
		
	}*/
	//int numconst = getNumConstraint();
	/*for(int i=0;i<numconst;i++)
	{
		//printf("zl[%d]:%lf\n",i,ZL[i]);
		//printf("zu[%d]:%lf\n",i,ZU[i]);
		z_l[i] = 1;
		z_u[i] = 1;
	}*/
}	

StepIImpl::~StepIImpl() {
	// TODO Auto-generated destructor stub
	delete[] irow;
	delete[] icol;
	delete[] ivalue;
}
/**
 * Evaluate the objective sum{i in RII} yi^(alpha)+1/(1-alpha)*sum{ j in 1 to |Vi|} pji^(alpha)<=Mialpha^{max}.
 * Input:
 *		x: Beamlet weight
 *
 *	@Author: Paras Babu Tiwari
 */

long double StepIImpl::evaluateObjective(const double *x) {
	long double sm = 0;

	assert(x != NULL);

	for (int i = 0; i < mtlb->getNumStepIIOar(); i++) {
		sm += x[YINDEX + i];
	}

	long double pjsm = 0;

	int k = 0;
	double **alpha = mtlb->getAlpha();
	//double *alpha = mtlb->getAlpha();
	//double alpha[1][1];
	//alpha[0][0] = 0.5;
	int *numvoxel = mtlb->getNumVoxel();
	for (int i = 0; i < mtlb->getNumStepIIOar(); i++) {
		int id = mtlb->getStepIId(i);
		int numVox = numvoxel[id];

		int index = 0;
		for (int j = 0; j < numVox; j++) {
			index = PINDEX + k;

			assert(index < getNumVar());
			pjsm += x[index];
			k++;
		}

		double factr = (1 - alpha[i][0] / 100) * numVox;
		double f = 1.0 / factr;
		pjsm = pjsm * f;

	}
	sm = sm + pjsm;
	return sm;

}
void StepIImpl::getBounds(double *x_l, double *x_u, double *g_l, double *g_u) {
	//printf("slip:%lf\n",Constant::SLIP);
	double slip = (1 + Constant::SLIP);
	printf("slip:%lf\n",slip);
	this->getBounds(x_l, x_u, g_l, g_u, slip);
}

/**
 * Get the upper and lower bounds of variables.
 * @Input Param: None
 * @Output Param:
 * 	x_l: lower bound for the variables.
 * 	x_u: Upper bound for the variables.
 * 	g_l: Lower bound for constraints.
 * 	g_u: Upper bound for constraints.
 */
void StepIImpl::getBounds(double *x_l, double *x_u, double *g_l, double *g_u,
		double slip,double s2) {

	//Set the upper and lower bound of infinity or max beamlet weight and 0 for variables.
	for (int i = 0; i < this->getNumVar(); i++) {
		x_l[i] = 0.0;
		if (i < NUMBEAMLET) {
			x_u[i] = mtlb->getUBWeight();
		} else {
			x_u[i] = 2e19;
		}
	}
	for (int i = 0; i < NUMZS; i++) {
		int index = ZINDEX+i;
		x_l[index] = -2e19;
		x_u[index] = 2e19;
	}
	int t_l = 0, t_u = 0;
	//Bound for Dj(w) - zji = 0.
	for (int i = 0; i < getNumFirstConst(); i++) {
		assert(t_l < getNumConstraint());
		assert(t_u < getNumConstraint());
		g_l[t_l] = 0;
		g_u[t_u] = 0;
		t_l++;
		t_u++;
	}
	//Bound for pjiα − zji + yiα>=0
	for (int i = 0; i < getNumSecondConst(); i++) {
		assert(t_l < getNumConstraint());
		assert(t_u < getNumConstraint());
		g_l[t_l] = 0;
		g_u[t_u] = 2e19;
		t_l++;

		t_u++;
	}
//	//Dj(w) ≤ Dimax.
	for (int i = 0; i < getNumThirdConst(); i++) {
		assert(t_l < getNumConstraint());
		assert(t_u < getNumConstraint());
		g_l[t_l] = -2e19;
		g_u[t_u] = mtlb->getMaxDose();
		t_l++;
		t_u++;
	}
//	// Dimin ≤ Dj(w) ≤ Dimax.
//	for (int i = 0; i < getNumFourthConst(); i++) {
//
//		assert(t_l < getNumConstraint());
//		assert(t_u < getNumConstraint());
//		double m;
//		m = step1->getMinDose();
//		//printf("Min dose:%lf\n", m);
//		g_l[t_l] = m *(1-s2);
//		printf("Min dose:%lf\n", g_l[t_l]);
//		m = step1->getMaxDose();
//		g_u[t_u] = m;
//		printf("Max dose:%lf\n", g_u[t_u]);
//		t_l++;
//		t_u++;
//
//	}
	int *numVoxel = mtlb->getNumVoxel();
	for(int i=0;i<mtlb->getNumTargets();i++)
	{
		int id = mtlb->getTargetsId(i);
		double mn,max;
		mn = step1->getMinDose(i);
		max = step1->getMaxDose(i);
		//printf("Min dose:%lf\n", mn);
		//printf("Max dose:%lf\n", max);
		for(int j=0;j<numVoxel[id];j++)
		{
			g_l[t_l] = mn *(1-s2);

			g_u[t_u] = max;

			t_l++;
			t_u++;
		}
	}
	//printf("Inside bounds\n");
	double *G = new double[mtlb->getNumTargets()];
	//printf("StepI :%x", step1->StepI());
	step1->evaluateG(G, step1->getW1());
//	// Gi(w) ≤ (1+s) Gi(wI).
	for (int i = 0; i < getNumFifthConst(); i++) {

		assert(t_l < getNumConstraint());
		assert(t_u < getNumConstraint());
		g_l[t_l] = -2e19;
		g_u[t_u] = G[i] * slip;
		//g_u[t_u]=0.0638017;
		printf("Upper bound:%lf\n",g_u[t_u]);
		t_l++;
		t_u++;
	}

	delete[] G;
	G = NULL;
	assert(t_l == getNumConstraint());
	assert(t_u == getNumConstraint());
}
/*
 * This function is a wrapper for the evalg.
 * It calls evalg() function adding total parameter.
 */
void StepIImpl::evalg(const double *x, double *g) {
	int total = 0;
	this->evalg(x, g, total);
	//printf("Total in evalg:%d\n",total);
	assert(total==getNumConstraint());

}
/*
 * Evaluate the constraint.
 * @Param:
 * 	 x: Beamlet weight.
 * 	 g: Output parameter, holds the constraints value at x.
 * 	@author:
 * 		Paras Babu Tiwari
 */
void StepIImpl::evalg(const double *x, double *g, int& total) {

	assert(x != NULL);
	assert(g != NULL);
	//Dj(w) - zji = 0.
	eval1stConst(x, total, g);
	eval2ndConst(x, total, g);
	eval3rdConst(x, total, g);
	eval4thConst(x, total, g);
	eval5thConst(x, total, g);

//	Matrices *D = new Matrices(mtlb->getTotalStepIIVoxels(), 1);
//	int *stepIIid = mtlb->getStepIId();
//
//	//Dj(w) - zji = 0.
//	mtlb->getDose(D, x, stepIIid, mtlb->getNumStepIIOar());
//
//	for (int i = 0; i < getNumFirstConst(); i++) {
//
//		assert(total < getNumFirstConst());
//		int index = ZINDEX + i;
//		assert(index < getNumVar());
//		assert(i < D->getRow());
//		g[total] = D->getAt(i, 0) - x[index];
//		total++;
//
//	}
//
//	delete D;
//
//	//pjiα − zji + yiα.
//	int numVox = 0;
//	for (int i = 0; i < mtlb->getNumStepIIOar(); i++) {
//
//		for (int j = 0; j < mtlb->getNumVoxels(stepIIid[i]); j++) {
//			for (int k = 0; k < mtlb->getNumVoxels(stepIIid[i]); k++) {
//				assert(total < getNumConstraint());
//				g[total] = x[PINDEX + j + i * numVox]
//						- x[ZINDEX + k + i * numVox] + x[YINDEX + i];
//				total++;
//			}
//		}
//		numVox = mtlb->getNumVoxels(stepIIid[i]);
//
//	}
//
//
//	//Dj(w) ≤ Dimax.
//	double *x1 = new double[NUMBEAMLET];
//	for (int i = 0; i < NUMBEAMLET; i++)
//		x1[i] = x[i];
//	Matrices *w1 = new Matrices(x1, NUMBEAMLET);
//
//	assert(step1 != NULL);
//	Matrices *DS = new Matrices(mtlb->getTTLStepIOar(), 1);
//	step1->getStepIDose(w1, DS);
//	for (int i = 0; i < DS->getRow(); i++) {
//		assert(i < getNumThirdConst());
//		assert(total < getNumConstraint());
//		g[total] = DS->getAt(i, 0);
//		total++;
//	}
//	delete DS;
//
//	//	// Dimin ≤ Dj(w) ≤ Dimax.
//	Matrices *DT = new Matrices(mtlb->getTTLTargetVoxel(), 1);
//	step1->getTargetDose(w1, DT);
//
//
//	for (int i = 0; i < getNumFourthConst(); i++) {
//		assert(total < getNumConstraint());
//		assert(i < DT->getRow());
//		g[total] = DT->getAt(i, 0);
//
//		total++;
//	}
//
//	delete DT;

	// Gi(w) ≤ (1+s) Gi(wI).
//	double *G = new double[mtlb->getNumTargets()];
//	step1->evaluateG(G);
//	for (int i = 0; i < mtlb->getNumTargets(); i++) {
//		assert(total <= getNumConstraint());
//		g[total] = G[i];
//		total++;
//	}
//	assert(total == getNumConstraint());
//	delete G;

//	for(int i=0;i<total;i++)
//	{
//		printf("g[%d]:%lf\n",g[i]);
//	}
//	printf("total:%d\n",total);

}

/*
 *	Evaluate the constraint Dj(w) - zji = 0.
 *	Input:
 *		x: Beamlet weight
 *		total: index in the constraint array.
 *	@Author: Paras Babu Tiwari
 */
void StepIImpl::eval1stConst(const double *x, int &total, double *g) {
	PrioritizedMain::Vector **D =
			new PrioritizedMain::Vector *[mtlb->getNumStepIIOar()];
	//Matrices *D = new Matrices(mtlb->getTotalStepIIVoxels(), 1);
	int *stepIIid = mtlb->getStepIId();

	//Dj(w) - zji = 0.
	mtlb->getDose(D, x, stepIIid, mtlb->getNumStepIIOar());

	//for (int i = 0; i < getNumFirstConst(); i++) {
	for (int i = 0; i < mtlb->getNumStepIIOar(); i++) {
		for (int j = 0; j < D[i]->size(); j++) {
			//assert(total < getNumFirstConst());
			int index = ZINDEX + i * D[i]->size() + j;
			//assert(index < getNumVar());
			//assert(i < D->getRow());
			//printf("index:%d\n",index);
			//g[total] = D->getAt(i, 0) - x[index];
			g[total] = D[i]->getAt(j) - x[index];
			total++;
		}

	}
	delete D;
}
/*
 *	Evaluate the constraint pjiα − zji + yiα.
 *	Input:
 *		x: Beamlet weight
 *		total: index in the constraint array.
 *	@Author: Paras Babu Tiwari
 */
void StepIImpl::eval2ndConst(const double *x, int &total, double *g) {
	int *stepIIid = mtlb->getStepIId();
	int numVox = 0;
	int *numvoxel = mtlb->getNumVoxel();
	for (int i = 0; i < mtlb->getNumStepIIOar(); i++) {
		int end = numvoxel[stepIIid[i]];
		for (int j = 0; j < end; j++) {
			assert(total < getNumConstraint());
			g[total] = x[PINDEX + j] - x[ZINDEX + j] + x[YINDEX + i];
			total++;
		}
	}

}

///*
// *	Evaluate the constraint pjiα − zji + yiα.
// *	Input:
// *		x: Beamlet weight
// *		total: index in the constraint array.
// *	@Author: Paras Babu Tiwari
// */
//void StepIImpl::eval2ndConst(const double *x, int &total, double *g) {
//	int *stepIIid = mtlb->getStepIId();
//	int numVox = 0;
//	for (int i = 0; i < mtlb->getNumStepIIOar(); i++) {
//
//		for (int j = 0; j < mtlb->getNumVoxels(stepIIid[i]); j++) {
//			for (int k = 0; k < mtlb->getNumVoxels(stepIIid[i]); k++) {
//				assert(total < getNumConstraint());
//				g[total] = x[PINDEX + j + i * numVox]
//						- x[ZINDEX + k + i * numVox] + x[YINDEX + i];
//				total++;
//			}
//		}
//		numVox = mtlb->getNumVoxels(stepIIid[i]);
//
//	}
//}

/*
 *	Evaluate the constraint Dj(w) <= Dimax.
 *	Input:
 *		x: Beamlet weight
 *		total: index in the constraint array.
 *	@Author: Paras Babu Tiwari
 */
void StepIImpl::eval3rdConst(const double *x, int &total, double *g) {
	double *x1 = new double[NUMBEAMLET];
	for (int i = 0; i < NUMBEAMLET; i++)
		x1[i] = x[i];
	PrioritizedMain::Vector *w1 = new PrioritizedMain::Vector(x1, NUMBEAMLET);

	assert(step1 != NULL);
	PrioritizedMain::Vector **DS =
			new PrioritizedMain::Vector *[mtlb->getNumStepIOar()];
	mtlb->getDose(DS, w1, mtlb->getStepIOrganId(), mtlb->getNumStepIOar());
	//Matrices *DS = new Matrices(mtlb->getTTLStepIOar(), 1);
	//step1->getStepIDose(w1, DS);
	for (int i = 0; i < mtlb->getNumStepIOar(); i++) {
		for (int j = 0; j < DS[i]->size(); j++) {
			//assert(i < getNumThirdConst());
			assert(total < getNumConstraint());
			g[total] = DS[i]->getAt(j);
			total++;
		}
	}
	delete w1;
	for (int i = 0; i < mtlb->getNumStepIOar(); i++) {
		delete DS[i];
	}
	delete[] DS;
}
/*
 *	Evaluate the constraint Dimin <= Dj(w) <= Dimax.
 *	Input:
 *		x: Beamlet weight
 *		total: index in the constraint array.
 *	@Author: Paras Babu Tiwari
 */
void StepIImpl::eval4thConst(const double *x, int &total, double *g) {
	//Matrices *DT = new Matrices(mtlb->getTTLTargetVoxel(), 1);
	PrioritizedMain::Vector *w1 = new PrioritizedMain::Vector(x, NUMBEAMLET);
	PrioritizedMain::Vector **DT =
			new PrioritizedMain::Vector *[mtlb->getNumTargets()];
	step1->getTargetDose(w1, DT);
	//printf("Target dose\n");
	//DT->printNonZeroEntry();
	for (int i = 0; i < mtlb->getNumTargets(); i++) {
		for (int j = 0; j < DT[i]->size(); j++) {
			assert(total < getNumConstraint());
			//assert(i < DT->getRow());
			g[total] = DT[i]->getAt(j);

			total++;
		}
	}
	delete w1;
	for (int i = 0; i < mtlb->getNumTargets(); i++) {
		delete DT[i];
	}
	delete[] DT;
}
/*
 *	Evaluate the constraint Gi(w) <= (1+s) Gi(wI).
 *	Input:
 *		x: Beamlet weight
 *		total: index in the constraint array.
 *	@Author: Paras Babu Tiwari
 */
void StepIImpl::eval5thConst(const double *x, int &total, double *g) {
	double *G = new double[mtlb->getNumTargets()];
	PrioritizedMain::Vector *w1 = new PrioritizedMain::Vector(x, NUMBEAMLET);
	step1->evaluateG(G, w1);
	for (int i = 0; i < mtlb->getNumTargets(); i++) {
		assert(total <= getNumConstraint());
		g[total] = G[i];
		total++;
	}
//assert(total == getNumConstraint());
	delete w1;
	delete[] G;
}
/*
 * Get the number of variables in the optimization.
 */
int StepIImpl::getNumVar() {

//We need to take into account of y's,p's and z's. I assumed that each organ has not more than one alpha value.So # y's equals to # of step2 organs.

//	printf("Influence Col:%d\n", mtlb->getColInfMat());
//	printf("# of Step II organ:%d\n", mtlb->getNumStepIIOar());
//	printf("Total # of voxels:%d\n", mtlb->getTotalStepIIVoxels());

	int n = NUMBEAMLET + NUMYS + NUMPS + NUMZS;
//printf("# of variables:%d\n",n);
	return n;

}
/*
 * @Input: None
 * @Return:
 * 	Number of non-zero entries in Jacobian matrix.
 * @Author:
 * 	Paras Babu Tiwari
 */
int StepIImpl::getNumJac() {

//pjiα − zji + yiα. Sum of total # of stepII voxels for p and z, and total number of step II organs to account for y.
//Dj(w) ≤ Dimax. Equals to total number of non-zero voxels in the step I organs.
// Gi(w) ≤ (1+s) Gi(wI). The # of these constraints equals to # of non-zero voxels in targets.
// Dimin ≤ Dj(w) ≤ Dimax. The # of these constraints equals to total # of non-zero voxels in target structures.
//	int n = mtlb->getNonzeroStepIIVoxels() + mtlb->getTotalStepIIVoxels()
//			+ mtlb->getTotalStepIIVoxels() + mtlb->getTotalStepIIVoxels()
//			+ mtlb->getNumStepIIOar() + mtlb->getNonzeroStepIOarVoxels()
//			+ mtlb->getNonZeroVoxelofTargets()
//			+ mtlb->getNonZeroVoxelofTargets();
//Dj(w) - zji = 0. The total number of non-zero entry from this constraint equals to the sum of
//# of non-zero entry in the influence matrix corresponding to the step II voxel, and # total of stepII organs voxel to account for z's.
	int numJacFrmFrstConst = getNumJacFromFirstConst();
	//printf("jac first const:%d\n", numJacFrmFrstConst);
//printf("Step II voxels:%d\n",mtlb->getTotalStepIIVoxels());
//Total # of second constraints
//int num2ndConst = mtlb->getTotalStepIIVoxels() * mtlb->getNumStepIIOar()
//		* mtlb->getNumStepIIOar() * mtlb->getTotalStepIIVoxels();
//Each constraint gives # of artifical variable entries in Jacobian matrix.
//	int numJacFrmSecondConst = num2ndConst * ARTIFICIALVAR;
	int numJacFrmSecondConst = getNumJacFromSecondConst();
	//printf("jac second const:%d\n", numJacFrmSecondConst);
	int numJacFromThirdConst = getNumJacFromThirdConst(); // mtlb->getNonzeroStepIOarVoxels();
	//printf("jac third const:%d\n", numJacFromThirdConst);
	int numJacFromFourthConst = getNumJacFromFourthConst(); // mtlb->getNonZeroVoxelofTargets();
	//printf("jac fourth const:%d\n", numJacFromFourthConst);
	int numJacFromFifthConst = getNumJacFromFifthConst();
	//printf("jac fifth const:%d\n", numJacFromFifthConst);
//printf("numJacFrmFrstConst:%d\n", numJacFrmFrstConst);
	//numJacFrmSecondConst = 0;
//	printf("numJacFrmFrstConst:%d\n",numJacFrmFrstConst);
//	printf("numJacFrmSecondConst:%d\n",numJacFrmSecondConst);
//	printf("numJacFromThirdConst:%d\n",numJacFromThirdConst);
//	printf("numJacFromFourthConst:%d\n",numJacFromFourthConst);
//	printf("numJacFromFifthConst:%d\n",numJacFromFifthConst);
	int n = numJacFrmFrstConst + numJacFrmSecondConst + numJacFromThirdConst
			+ numJacFromFourthConst + numJacFromFifthConst;

	return n;
}
/*
 * Get the number of constraints.
 *
 */
int StepIImpl::getNumConstraint() {

// Gi(w) ≤ (1+s) Gi(wI). The # of these constraints equals to # of targets.
// Dimin ≤ Dj(w) ≤ Dimax. The # of these constraints equals to total # of voxels in target structures.

//	int numFirst = mtlb->getTotalStepIIVoxels();
//
//	int numSecond = mtlb->getTotalStepIIVoxels() * mtlb->getNumStepIIOar()
//			* mtlb->getNumStepIOar() * mtlb->getTotalStepIIVoxels();
//	+ mtlb->getTTLStepIOar() + mtlb->getNumTargets()
//	+ mtlb->getTTLTargetRow();
	int m1, m2, m3, m4, m5;
	m1 = getNumFirstConst();
	m2 = getNumSecondConst();
	m3 = getNumThirdConst();
	m4 = getNumFourthConst();
	m5 = getNumFifthConst();

	int m = m1 + m2 + m3 + m4 + m5;
//printf("1st const:%d 2nd const:%d 3rd const:%d 4th const:%d\n",getNumFirstConst(),getNumSecondConst(),getNumThirdConst(),getNumFourthConst());
	return m;

}
/**
 * Dj(w) - zji = 0. There will be getTotalStepIIVoxels()(Total # of voxels in step 2 organs) number of constraints of this type.
 *	@return the # of Dj(w) - zji = 0 constraints
 *	@Author: Paras Babu Tiwari
 */
int StepIImpl::getNumFirstConst() {

	return mtlb->getTotalStepIIVoxels();

}
/**
 *	pjiα − zji + yiα. There will be #of voxel * # of StepII organ * # of stepII organ * # of voxel constraints. For each p there is # of voxel
 choice for z, that gives the first term. There is # of stepII organ choice for y for a give p and z, which gives the second term.
 For each voxel of z, there is a choice of # of organs, which gives the third term. Finally, p's index can go for # of voxel time.
 */
int StepIImpl::getNumSecondConst() {
//	int numSecond = mtlb->getTotalStepIIVoxels() * mtlb->getNumStepIIOar()
//			* mtlb->getNumStepIIOar() * mtlb->getTotalStepIIVoxels();
//
	//return 0;
//	int numSecond = 0;
//	int numVoxel = 0;
//	for (int i = 0; i < mtlb->getNumStepIIOar(); i++) {
//
//		int id = mtlb->getStepIId(i);
//		numVoxel = numvoxel[id];
//		numSecond += numVoxel * numVoxel;
//
//	}
//	return numSecond;
	return mtlb->getTotalStepIIVoxels();
}
/*
 * 	Dj(w) ≤ Dimax. Equals to total number of voxels in the step I organs.
 * 	@return: # of voxels in the step I organ.
 * 	@Author: Paras Babu Tiwari
 */
int StepIImpl::getNumThirdConst() {

	//return 0;
	return mtlb->getTTLStepIOar();
}
/**
 * Dimin ≤ Dj(w) ≤ Dimax. The # of these constraints equals to total # of target voxels.
 * @return: Total # of constraints
 *
 * @Author: Paras Babu Tiwari
 *
 */
int StepIImpl::getNumFourthConst() {
	//return 0;

	return mtlb->getTTLTargetVoxel();

}
/**
 * Gi(w) ≤ (1+s) Gi(wI). The # of these constraints equals to # of targets.
 */
int StepIImpl::getNumFifthConst() {
	//return 0;

	return mtlb->getNumTargets();

}
/*
 * Dj(w) - zji = 0.The getNonzeroStepIIVoxels() gives the total # of voxel that has non-zero entry
 * in the influence matrix which accounts for the contribution from Dj(w). And total # of voxels is contributed by the z term.
 * @return: Total # of entry in the Jacobian matrix.
 * @Author Paras Babu Tiwari
 */
int StepIImpl::getNumJacFromFirstConst() {

	return mtlb->getNonzeroStepIIVoxels() + mtlb->getTotalStepIIVoxels();
}
/**
 * pjiα − zji + yiα.
 */
int StepIImpl::getNumJacFromSecondConst() {

	int numSecond = 0;
//	int numVoxel = 0;
//	for (int i = 0; i < mtlb->getNumStepIIOar(); i++) {
//
//		int id = mtlb->getStepIId(i);
//		numVoxel = numvoxel[id];
//		numSecond += numVoxel * numVoxel * ARTIFICIALVAR;
//
//	}
	numSecond = mtlb->getTotalStepIIVoxels() * ARTIFICIALVAR;
	return numSecond;

}
/**
 * Dj(w) ≤ Dimax.The number of Jacobian entries equals to the # of voxels that has non-zero entry in the influence matrix.
 * @return: # of entries in the Jacobian matrix.
 * @Author: Paras Babu Tiwari
 */
int StepIImpl::getNumJacFromThirdConst() {

	//return 0;
	return mtlb->getNonzeroStepIOarVoxels();
}
/**
 * Dimin ≤ Dj(w) ≤ Dimax. The # of these constraints equals to total # of target voxels that has non-zero entry in the influence matrix ..
 * @return: Total # of entries in Jacobian matrix due to this constraint.
 *
 * @Author: Paras Babu Tiwari
 *
 */
int StepIImpl::getNumJacFromFourthConst() {

	//return 0;
	return mtlb->getNonZeroVoxelofTargets();

}
/**
 * Dimin ≤ Dj(w) ≤ Dimax. The # of these constraints equals to total # of target voxels that has non-zero entry in the influence matrix ..
 * @return: Total # of entries in Jacobian matrix due to this constraint.
 *
 * @Author: Paras Babu Tiwari
 *
 */
int StepIImpl::getNumJacFromFifthConst() {


	if(numFifth>0)
		return numFifth;
	//return 0;
	//return mtlb->getNonZeroVoxelofTargets();
	//return NUMBEAMLET;
	//printf("Inside from fith constraint\n");
	int count = 0;
	int *numvoxel = mtlb->getNumVoxel();
	for (int i = 0; i < mtlb->getNumTargets(); i++) {
		int id = mtlb->getTargetsId(i);
		int n = numvoxel[id];
		double *row = mtlb->getVoxels()[id];
		for (int j = 0; j < NUMBEAMLET; j++) {
			for (int k = 0; k < n; k++) {
				int r = row[k];
				double val = 0;

				double coeff = mtlb->influence->getAt(r, j);

				if (coeff != 0) {
					count++;
					break;
				}
			}
		}
	}
	numFifth = count;
	//printf("done fifth\n");
	return count;

}
/**
 * Evaluate the gradient of function.
 * The function F(x) = yi(alpha) + 1/(1-alpha)*|vi| sum{j=1 to |vi|} pji(alpha).
 * The gradient of a function is given by: \/f(x) = (df/dx1, df/dx2,df/dx3...df/dxn df/dy1 df/dy2..df/dp1 df/dp2...df/z1 df/dz2..)
 * The gradient w.r.to x's are zero as there is no x term in the objective function.
 * The gradient w.r.to y's equals to one as y's terms are linear.
 * The gradient w.r.to p's equals to 1/(1-alpha)*|vi| because p's are also linear.
 * The gradient w.r.to z's are zero as there is no z term in the objective function.
 */
void StepIImpl::evaluategradf(const double *x, double *g) {

	double factor;
	assert(g != NULL);
	for (int i = 0; i < getNumVar(); i++)
		g[i] = 0;
//The number of beamlets equals to the column of the influence matrix. So the YINDEX must be passed the # of column of Influence Matrices.
	int YINDEX = mtlb->getColInfMat();
//printf("YINDEX:%d\n", YINDEX);
	int index = 0;
	for (int i = 0; i < mtlb->getNumStepIIOar(); i++) {
		index = YINDEX + i;
		assert(index < getNumVar());
		g[index] = 1;
	}
//p's are after beamlets and y's in the array.

//printf("pindex:%d\n", PINDEX);
	int *numvoxel = mtlb->getNumVoxel();
	double **alpha = mtlb->getAlpha();
	//double *alpha = mtlb->getAlpha();

	for (int i = 0; i < mtlb->getNumStepIIOar(); i++) {
		int id = mtlb->getStepIId(i);
		int numVox = numvoxel[id];
		factor = (1 - alpha[i][0] / 100) * numVox;
		factor = 1.0 / factor;
		for (int j = 0; j < numVox; j++) {
			index = PINDEX + j;
			assert(index < getNumVar());
			g[index] = factor;
		}
		//printf("pindex+i:%d\n", PINDEX + i);
	}
}
void StepIImpl::getJacobian(int *iRow, int *iCol, double *value,
		const double *x) {
	int total = 0;
	int l = 0;
	this->getJacobian(iRow, iCol, value, x, total, l);
	//printf("total in jacobian:%d\n",total);
	assert(total == getNumJac());

}
/*
 * Return either the sparsity structure of the Jacobian of the constraints, or the values for the Jacobian of the constraints at the point x.
 * The detail about Jacobian matrix can be found at http://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant
 * @Input:
 * 	iRow: Output Parameter, holds the row # of non-zero entry in the Jacobian matrix
 * 	iCol: Output Parameter, holds the column # of non-zero entry in the Jacobian matrix.
 * 	value: Output Parameter, holds the value of non-zero entry of the Jacobian matrix. The value[i] corresponds to iRow[i] and iCol[i] entry in the matrix.
 * 	x: Input Parameter: The solution vector found during the optimization.
 *  @Author: Paras Babu Tiwari
 *
 */
void StepIImpl::getJacobian(int *iRow, int *iCol, double *value,
		const double *x, int& total, int& l) {

	time_t start_t, end_t;
	start_t = time(NULL);

	int noOfConst = 0;
	if (value == NULL) {

		//Dj(w) - zji = 0.
		handleFirstConst(iRow, iCol, value, x, l, total);
		assert(l == getNumFirstConst());
		assert(total == getNumJacFromFirstConst());

//		//pjiα − zji + yiα>=0.
		handleSecondConst(iRow, iCol, value, x, l, total);
		noOfConst = getNumFirstConst() + getNumSecondConst();
		assert(l == noOfConst);
		int noOfJac = 0;
		noOfJac = getNumJacFromFirstConst() + getNumJacFromSecondConst();
		assert(total == noOfJac);

//		//Dj(w) <= Dimax.
		handleThirdConst(iRow, iCol, value, x, l, total);
		noOfJac = getNumJacFromFirstConst() + getNumJacFromSecondConst()
				+ getNumJacFromThirdConst();
		assert(total == noOfJac);
		noOfConst = getNumFirstConst() + getNumSecondConst()
				+ getNumThirdConst();
		assert(l == noOfConst);

		// Dimin ≤ Dj(w) ≤ Dimax.
		handleFourthConst(iRow, iCol, value, x, l, total);
		noOfConst = getNumFirstConst() + getNumSecondConst()
				+ getNumThirdConst() + getNumFourthConst();
		assert(l == noOfConst);
		noOfJac = getNumJacFromFirstConst() + getNumJacFromSecondConst()
				+ getNumJacFromThirdConst() + getNumJacFromFourthConst();
		assert(total == noOfJac);

		// Gi(w) ≤ (1+s) Gi(wI).
		handleFifthConst(iRow, iCol, value, x, l, total);
		assert(l == getNumConstraint());
	} //End of if
	else {

		//Dj(w) - zji = 0.
		handleFirstConst(iRow, iCol, value, x, l, total);
//
//			//pjiα − zji + yiα>=0.
		handleSecondConst(iRow, iCol, value, x, l, total);
//
//			//Dj(w) ≤ Dimax.
		handleThirdConst(iRow, iCol, value, x, l, total);
//
//			// Dimin ≤ Dj(w) ≤ Dimax.
		handleFourthConst(iRow, iCol, value, x, l, total);

		// Gi(w) ≤ (1+s) Gi(wI).
		handleFifthConst(iRow, iCol, value, x, l, total);

	}

	end_t = time(NULL);
	double diff = difftime(end_t, start_t);
	//printf("Elapsed time in Jacobian :%lf seconds\n", diff);

}
void StepIImpl::handleFourthConst(int *iRow, int *iCol, double *value,
		const double *x, int &l, int &total) {

	//The Jacobian matrix entries are zero for y,p and z variables as the constraint doesn't has these variables.
	// Jac(i) = [ai1 ai2 ai3..ain...0...0]
	double val = 0;
	struct timeval start, end;
	long mtime, seconds, useconds;
	int *numvoxel = mtlb->getNumVoxel();	
	for (int i = 0; i < numTargets; i++) 
	{
		int oldtotal = total;
		int id = mtlb->getTargetsId(i);
		int n = numvoxel[id];
		double *row = mtlb->getVoxels()[id];

		if (fjacsz[i] == 0) {
			Matrices *res;			
			res = mtlb->influence->getSubMatrix(id, mtlb->voxels[id],
					numvoxel[mtlb->getTargetsId(i)]);
			std::vector<double> nzval = res->getNzElements();			
			fjacsz[i] = nzval.size();
			
			fourthnumrow[i] = res->getRow();
			std::vector<long int> rows = res->getNzRows();
			std::vector<long int> cols = res->getNzCols();
			fourthconstjac[i] = new double[fjacsz[i]];			
			memcpy(fourthconstjac[i], &nzval[0], sizeof(double) * fjacsz[i]);
			delete res;
		}		
		if (value != NULL) {
			memcpy(value + total, fourthconstjac[i], sizeof(double) * fjacsz[i]);
			total = total + fjacsz[i];
			l = l + fourthnumrow[i];
		} 
		else 
		{
			Matrices *res;			
			res = mtlb->influence->getSubMatrix(id, mtlb->voxels[id],
					numvoxel[mtlb->getTargetsId(i)]);
			std::vector<long int> rows = res->getNzRows();
			std::vector<long int> cols = res->getNzCols();
			for (int k = 0; k < fjacsz[i]; k++) 
			{
				iRow[total] = l + rows[k];
				iCol[total] = cols[k];
				total++;
			}
			l = l + fourthnumrow[i];
			delete res;			
			

		}

	} //End of for
	//printf("l at the end of foruth :%d\n", l);

} 
/*
void StepIImpl::handleFourthConst(int *iRow, int *iCol, double *value,
		const double *x, int &l, int &total) {

//The Jacobian matrix entries are zero for y,p and z variables as the constraint doesn't has these variables.
// Jac(i) = [ai1 ai2 ai3..ain...0...0]
	double val = 0;
//printf("l in fourth :%d\n", l);
	int *numvoxel = mtlb->getNumVoxel();
	for (int i = 0; i < mtlb->getNumTargets(); i++) {
		int id = mtlb->getTargetsId(i);
		int n = numvoxel[id];
		double *row = mtlb->getVoxels()[id];

		for (int k = 0; k < n; k++) {
			//printf("k:%d\n",k);
			int r = row[k];
			for (int j = 0; j < NUMBEAMLET; j++) {
				val = mtlb->influence->getAt(r, j);
				//Fill in the matrix structure as value is null.
				if (val != 0 && value == NULL) {
					assert(total < getNumJac());

					iRow[total] = l; //r;
					iCol[total] = j; //j;
//					irow[total] = l;
//					icol[total] = j;
							//printf("l:%d j:%d\n",l,j);
					total++;
				} //End of if
				else if (val != 0 && value != NULL) {
					assert(total < getNumJac());
					//Set the value.
					value[total] = val;
					//ivalue[total] = val;
					total++;
				}
			} //End of for
			l++;
		} //End of for
	} //End of for
// printf("l at the end of foruth :%d\n", l);

} //End of function
*/
/**
 *	 Set either row,column indices or the value array for the constraint  Gi(w) ≤ (1+s) Gi(wI).
 *	 @Input:
 * 		iRow: Output Parameter, holds the row # of non-zero entry in the Jacobian matrix
 * 		iCol: Output Parameter, holds the column # of non-zero entry in the Jacobian matrix.
 * 		value: Output Parameter, holds the value of non-zero entry of the Jacobian matrix. The value[i] corresponds to iRow[i] and iCol[i] entry in the matrix.
 * 		x: Input Parameter: The solution vector found during the optimization.
 * 	@Author: Paras Babu Tiwari
 */
/**
 *	 Set either row,column indices or the value array for the constraint Dimin ≤ Dj(w) ≤ Dimax.
 *	 @Input:
 * 		iRow: Output Parameter, holds the row # of non-zero entry in the Jacobian matrix
 * 		iCol: Output Parameter, holds the column # of non-zero entry in the Jacobian matrix.
 * 		value: Output Parameter, holds the value of non-zero entry of the Jacobian matrix. The value[i] corresponds to iRow[i] and iCol[i] entry in the matrix.
 * 		x: Input Parameter: The solution vector found during the optimization.
 * 	@Author: Paras Babu Tiwari
 */
 void StepIImpl::handleFifthConst(int *iRow, int *iCol, double *value,
		const double *x, int &l, int &total) {
	double coeff = 0;
	long double val = 0;
	int n = 0;
	int numTarget = numTargets;
	PrioritizedMain::Vector **D = NULL;
	Matrices *res;
	int *numvoxel = mtlb->getNumVoxel();
	struct timeval start, end;
	long mtime, seconds, useconds;	
	PrioritizedMain::Vector **jac = new PrioritizedMain::Vector *[numTarget];	
	if (value != NULL) {
		D = new PrioritizedMain::Vector*[numTarget];
		step1->getTargetDose(x, D);	
		step1->subPrescDose(D);		
		step1->divideByVoxel(D);
	}	
	for (int i = 0; i < numTargets; i++) {
		if (value != NULL) {
			D[i] = (*D[i]) * 2;			
			jac[i] = targetInf[i]->matrixTransMultiply(D[i]);			
			int sz = jac[i]->size();			
			for (int j = 0; j < sz; j++) {
				double val = jac[i]->getAt(j);				
				if (val != 0) {
					value[total] = val;
					total++;
				}
			}			
		} else {			
			for (int k = 0; k < targetInf[i]->getCol(); k++) {
				for (int j = 0; j < targetInf[i]->getRow(); j++) {
					double val = targetInf[i]->getAt(j, k);
					if (val != 0) {
						iRow[total] = l;
						iCol[total] = k;
						total++;
						break;
					}
				}
			}			
		}		
		l++;
	}	
	if (D != NULL) {
		for (int i = 0; i < numTargets; i++) {
			delete D[i];
		}
		delete[] D;
	}
	

}
 /*
void StepIImpl::handleFifthConst(int *iRow, int *iCol, double *value,
		const double *x, int &l, int &total) {

//The Jacobian matrix entries are zero for y,p and z variables as the constraint doesn't has these variables.
// Jac(i) = 2/|Vt| *(Dj(w)-Dipre)*[ai1 ai2 ai3..ain...0...0]

	double coeff = 0;
	long double val = 0;
	int n = 0;
	PrioritizedMain::Vector **D = NULL;
	if (value != NULL) {
		PrioritizedMain::Vector *w = new PrioritizedMain::Vector(x, NUMBEAMLET);
		D = new PrioritizedMain::Vector*[mtlb->getNumTargets()];
		step1->getTargetDose(w, D);
		step1->subPrescDose(D);
		//printf("D\n");
		//D->printNonZeroEntry();
		step1->divideByVoxel(D);
		delete w;
//		printf("D\n");
//		D->printNonZeroEntry();
	}

	int *numvoxel = mtlb->getNumVoxel();
	for (int i = 0; i < mtlb->getNumTargets(); i++) {
		int id = mtlb->getTargetsId(i);
		n = numvoxel[id];
		double *row = mtlb->getVoxels()[id];
		for (int j = 0; j < NUMBEAMLET; j++) {
			val = 0;
			for (int k = 0; k < n; k++) {
				int r = row[k];


				coeff = mtlb->influence->getAt(r, j);

				if (coeff != 0 && value == NULL) {

					//Fill in the matrix structure as value is null.

					iRow[total] = l;
					iCol[total] = j;
					total++;
					break;
				} else if (coeff != 0 && value != NULL) {

//					printf("t:%d\n",t);
//					printf("total:%d\n",total);
					//val = 2 * D->multiplyAt(t, 0, coeff);
					double t = 0;
					t = D[i]->getAt(k) * 2;
//					if (val == -1)
//						val = t * coeff;
//					else
					val += t * coeff;
					//val=0;
					//Set the value.
					//value[total] = val;
					//total++;
				} //End of else if
			} //End of for
			if (value != NULL && val != 0) {
				value[total] = val;
				total++;
			}

		} //End of for
		l++;

	} //End of for
	if (D != NULL) {
		for (int i = 0; i < mtlb->getNumTargets(); i++) {
			delete D[i];

		}
		delete[] D;
	}
// printf("l at the end of fifth :%d\n", l);
}*/
/**
 *	 Set either row,column indices or the value array for the constraint Dj(w) ≤ Dimax.
 *	 @Input:
 * 		iRow: Output Parameter, holds the row # of non-zero entry in the Jacobian matrix
 * 		iCol: Output Parameter, holds the column # of non-zero entry in the Jacobian matrix.
 * 		value: Output Parameter, holds the value of non-zero entry of the Jacobian matrix. The value[i] corresponds to iRow[i] and iCol[i] entry in the matrix.
 * 		x: Input Parameter: The solution vector found during the optimization.
 * 	@Author: Paras Babu Tiwari
 */
void StepIImpl::handleThirdConst(int *iRow, int *iCol, double *value,
		const double *x, int &l, int &total) {

	double val;
//printf("********* Inside third constraint **************\n");
//Dj(w) ≤ Dimax. Jac(i) = [ai1 ai2 .. ain..0..0 for all other variables.
	int *numvoxel = mtlb->getNumVoxel();
	for (int i = 0; i < mtlb->getNumStepIOar(); i++) {
		int id = mtlb->getStepIOrganId(i);
		//printf("id:%d\n",id);
		//Get the # of Voxels in the organ.
		int n = numvoxel[id];
		//printf("n:%d\n",n);
		//row[i] corresponds to ith row in influence matrix for the ith voxel of the organ having id 'id'.
		double *row = mtlb->getVoxels()[id];

		for (int k = 0; k < n; k++) {
			int r = row[k];

			for (int j = 0; j < NUMBEAMLET; j++) {
				//printf("row:%d col:%d:lf\n",r,j,val);
				val = mtlb->influence->getAt(r, j);
				if (val != 0) {
					//Fill in the matrix structure as value is null.
					if (value == NULL) {
						//printf("val:%lf\n",val);

						iRow[total] = l;
						iCol[total] = j;
//						irow[total] = l;
//						icol[total] = j;
						total++;
						//printf("l:%d j:%d\n",l,j);
					} else {
						//Set the value.
						value[total] = val;
						//ivalue[total] = val;
						total++;
					}

				} //End of if
			}
			l++;
		} //End of for

	} //End of for
//printf("l at the end of third :%d\n", l);
}
/**
 *	 Set either row,column indices or the value array for the constraint pjiα − zji + yiα>=0.
 *	 @Input:
 * 		iRow: Output Parameter, holds the row # of non-zero entry in the Jacobian matrix
 * 		iCol: Output Parameter, holds the column # of non-zero entry in the Jacobian matrix.
 * 		value: Output Parameter, holds the value of non-zero entry of the Jacobian matrix. The value[i] corresponds to iRow[i] and iCol[i] entry in the matrix.
 * 		x: Input Parameter: The solution vector found during the optimization.
 * 	@Author: Paras Babu Tiwari
 */
void StepIImpl::handleSecondConst(int *iRow, int *iCol, double *value,
		const double *x, int &l, int &total) {

//printf("Handle Second Constraint\n");
	int *numvoxel = mtlb->getNumVoxel();
//pjiα − zji + yiα>=0. [0..for w's terms..1 for y and p's terms, and -1 for z's terms.
	for (int i = 0; i < mtlb->getNumStepIIOar(); i++) {

		int id = mtlb->getStepIId(i);
		int numVoxel = numvoxel[id];
		for (int j = 0; j < numVoxel; j++) {

			if (value == NULL) {
				//Fill in the matrix structure as value is null.

				assert(total < getNumJac());

				//Entry for y
				iRow[total] = l;
				iCol[total] = YINDEX + i;
				//					irow[total] = l;
				//					icol[total] = j;
				total++;
				assert(total < getNumJac());
				//Entry for p

				iRow[total] = l;
				iCol[total] = PINDEX + j;
				//					irow[total] = l;
				//					icol[total] = j;
				total++;
				assert(total < getNumJac());

				//Entry for z
				iRow[total] = l;
				iCol[total] = ZINDEX + j;
				//					irow[total] = l;
				//					icol[total] = j;
				total++;
			} else if (value != NULL) {
				//Set the value.
				//set -1 for z.

				assert(total < getNumJac());
				//set 1 for y and p's
				value[total] = 1;
				//ivalue[total] = 1;
				total++;
				assert(total < getNumJac());
				value[total] = 1;
				//ivalue[total] = 1;
				total++;
				assert(total < getNumJac());
				value[total] = -1;
				//ivalue[total] = -1;
				total++;

			}

			l++;
		} //End of for

	} //End of for
// printf("l at the end of second :%d\n", l);
}
///**
// *	 Set either row,column indices or the value array for the constraint pjiα − zji + yiα>=0.
// *	 @Input:
// * 		iRow: Output Parameter, holds the row # of non-zero entry in the Jacobian matrix
// * 		iCol: Output Parameter, holds the column # of non-zero entry in the Jacobian matrix.
// * 		value: Output Parameter, holds the value of non-zero entry of the Jacobian matrix. The value[i] corresponds to iRow[i] and iCol[i] entry in the matrix.
// * 		x: Input Parameter: The solution vector found during the optimization.
// * 	@Author: Paras Babu Tiwari
// */
//void StepIImpl::handleSecondConst(int *iRow, int *iCol, double *value,
//		const double *x, int &l, int &total) {
//
//	//printf("Handle Second Constraint\n");
//
//	//pjiα − zji + yiα>=0. [0..for w's terms..1 for y and p's terms, and -1 for z's terms.
//	for (int i = 0; i < mtlb->getNumStepIIOar(); i++) {
//
//		int id = mtlb->getStepIId(i);
//		int numVoxel = numvoxel[id];
//		for (int j = 0; j < numVoxel; j++) {
//			for (int k = 0; k < numVoxel; k++) {
//
//				if (value == NULL) {
//					//Fill in the matrix structure as value is null.
//
//					assert(total < getNumJac());
//
//					//Entry for y
//					iRow[total] = l;
//					iCol[total] = YINDEX + i;
////					irow[total] = l;
////					icol[total] = j;
//					total++;
//					assert(total < getNumJac());
//					//Entry for p
//
//					iRow[total] = l;
//					iCol[total] = PINDEX + j;
////					irow[total] = l;
////					icol[total] = j;
//					total++;
//					assert(total < getNumJac());
//
//					//Entry for z
//					iRow[total] = l;
//					iCol[total] = ZINDEX + k;
////					irow[total] = l;
////					icol[total] = j;
//					total++;
//				} else if (value != NULL) {
//					//Set the value.
//					//set -1 for z.
//
//					assert(total < getNumJac());
//					//set 1 for y and p's
//					value[total] = 1;
//					//ivalue[total] = 1;
//					total++;
//					assert(total < getNumJac());
//					value[total] = 1;
//					//ivalue[total] = 1;
//					total++;
//					assert(total < getNumJac());
//					value[total] = -1;
//					//ivalue[total] = -1;
//					total++;
//
//				}
//
//				l++;
//			} //End of for
//
//		} //End of for
//	} //End of for
//	  // printf("l at the end of second :%d\n", l);
//}
/*
 * Set either row,column indices or the value array for the constraint Dj(w) - zji = 0.
 * @Input:
 * 	iRow: Output Parameter, holds the row # of non-zero entry in the Jacobian matrix
 * 	iCol: Output Parameter, holds the column # of non-zero entry in the Jacobian matrix.
 * 	value: Output Parameter, holds the value of non-zero entry of the Jacobian matrix. The value[i] corresponds to iRow[i] and iCol[i] entry in the matrix.
 * 	x: Input Parameter: The solution vector found during the optimization.
 * 	@Author: Paras Babu Tiwari
 *
 */
void StepIImpl::handleFirstConst(int *iRow, int *iCol, double *value,
		const double *x, int &l, int &total) {

//printf("Handle first constraint\n");
//printf("numPBs:%d\n",numPBs);
//	int n = getNumJac() - getNumJacFromFifthConst();
	double val = 0;
//	if (irow == NULL) {
//		irow = new int[n];
//	}
//	if (icol == NULL) {
//		icol = new int[n];
//	}
//	if (value != NULL && ivalue == NULL) {
//		ivalue = new double[n];
//	}
	int *numvoxel = mtlb->getNumVoxel();
//Dj(w) - zji = 0.Jac(i) = [ ai1 ai2 ai3... 0 .. upto y and p's terms, and -1 for z's terms].
	for (int i = 0; i < mtlb->getNumStepIIOar(); i++) {
		int id = mtlb->getStepIId(i);
		//Get the # of Voxels in the organ.
		int numVoxel = numvoxel[id];
		//printf("Numvoxel:%d\n",numVoxel);

		//row[i] corresponds to ith row in influence matrix for the ith voxel of the organ having id 'id'.
		double *row = mtlb->getVoxels()[id];

		for (int k = 0; k < numVoxel; k++) {

			int r = row[k];

//			printf("numPBs:%d\n",numPBs);
			//Substract the # of y's and p's as there entry is zero in Jacobian matrix
//			int totItr = this->getNumVar()
//					- (mtlb->getNumStepIIOar() + mtlb->getTotalStepIIVoxels());

			//There are numPBs beamlets, and one more to account for z.
			//for (int j = 0; j <= NUMBEAMLET; j++) {
			for (int j = 0; j <= NUMBEAMLET; j++) {
				//set w's terms
				if (j < NUMBEAMLET) {
					val = mtlb->influence->getAt(r, j);

					if (val != 0) {

						//Fill in the matrix structure in value is null.
						if (value == NULL) {

							assert(total < getNumJac());
							//assert(l < getNumFirstConst());
							iRow[total] = l;
							iCol[total] = j;
//							irow[total] = l;
//							icol[total] = j;
							total++;
						} else {
							assert(total < getNumJac());
							//Set the value.
							value[total] = val;
							//ivalue[total] = val;
							total++;
						}

					} //End of if
				} //End of if
				  //set z terms
				else {
					//Fill in the matrix structure if value is null.
					if (value == NULL) {
//						int k = j + mtlb->getNumStepIIOar()
//								+ mtlb->getTotalStepIIVoxels();
//						printf("k:%d\n",k);
						//printf("row:%d col:%d\n",l,k+l);
						assert(total < getNumJac());
						//assert(l < getNumFirstConst());
						iRow[total] = l;
						iCol[total] = ZINDEX + l;
//						irow[total] = l;
//						icol[total] = ZINDEX + l;
						total++;
					} else {
						assert(total < getNumJac());
						//Set the value.
						value[total] = -1;
						//ivalue[total] = -1;
						total++;
					}
				} //End of else

			} //End of for
			l++;
		} //End of for

	} //End of for
//printf("l at the end of first :%d\n", l);

}
/*
 * Get number of element in the hessian matrix.
 * @return number of non-zero element in the hessian Matrices
 */
int StepIImpl::getNumHess() {

	//return 0;
	if (mtlb->getHMatrix() == NULL)
		return 0;
//return mtlb->getHMatrix()->getNzEntry() - mtlb->getNumTargets();
#ifdef dummy
	return mtlb->getHMatrix()->getNzEntry();
#else
	return mtlb->getHMatrix()->getNzEntry() - mtlb->getNumTargets();
#endif

}
/*
 * Print the dose
 */
void StepIImpl::printDose() {

	printDose(w2);
}
/*
 * Print the dose obtained from the optimization
 */
void StepIImpl::printDose(PrioritizedMain::Vector *w) {
	int *stepIIid = mtlb->getStepIId();
	PrioritizedMain::Vector **D =
			new PrioritizedMain::Vector*[mtlb->getNumStepIIOar()];

	mtlb->getDose(D, w, stepIIid, mtlb->getNumStepIIOar());

	int start = 0;
	int end = 0;
	int id = 0, i = 0;
	double v;
	int *numvoxel = mtlb->getNumVoxel();
	for (i = 0; i < mtlb->getNumStepIIOar(); i++) {
		id = mtlb->getStepIId(i);
		end = start + numvoxel[id];
		v = D[i]->getMin();
		printf("Minimum Dose:%lf\n", v);
		v = D[i]->getMax();
		printf("Maximum Dose:%lf\n", v);
		v = D[i]->average();
		printf("Average Dose:%lf\n", v);
		start = end;
	}
	for (i = 0; i < mtlb->getNumStepIIOar(); i++) {
		delete D[i];
	}
	delete[] D;
//	char *fileName;
//	fileName = new char[FILESIZE];
//	sprintf(fileName, "w2%d", Constant::CASEID);
//	fileName[FILESIZE - 1] = '\0';
//	FileManager *f = new FileManager(fileName);
//	double *x;
//	x = new double[NUMBEAMLET];
//	w2->getDouble(x);
//	f->writeFile(x, NUMBEAMLET);
//	f->close();
}
bool StepIImpl::getHessian(int *iRow, int *iCol, double *value,
		double obj_factor, const double *lambda) {
	int total = 0;
	getHessian(iRow, iCol, value, obj_factor, lambda, total, false);
	return true;

}
/**
 * Get the Hessian. We have stored the hessian in matlab index format, so need to substract one.
 * Hessian is A'*A
 */
bool StepIImpl::getHessian(int *iRow, int *iCol, double *value,
		double obj_factor, const double *lambda, int &total, bool step4) {

	time_t start_t, end_t;
	time_t hstime, hendtime;
	start_t = time(NULL);

	int rowindex;
	int colindex;
	double val = 0;
	double totalhash = 0;
	if (mtlb->getHMatrix() == NULL)
		return false;
//printf("NUMBEAMLET:%d\n",NUMBEAMLET);
	if (value == NULL) {

		for (int i = 0; i < mtlb->getHMatrix()->getNzEntry(); i++) {
			//for (int i = 0; i < 280593; i++) {
			if (!step4)
				assert(total <= getNumHess());
			rowindex = (int) mtlb->getHMatrix()->getRowNumber(i); //mtlb->getHMatrix()->rowIndex[i];
			colindex = (int) mtlb->getHMatrix()->getColNumber(i); //mtlb->getHMatrix()->colIndex[i];

			if (rowindex < NUMBEAMLET) {
				//iRow[total] = rowindex - 1;
				//iCol[total] = colindex - 1;
				//printf("rowIndex:%d colIndex:%d\n",rowindex,colindex);
				//Step4 also call this method to fill the hessian structure.
				//We already specified the row number and colnumber for the diagonal element in stepIV class, so no need to do again here.
				if (step4) {
					if (rowindex == colindex) {
						continue;
					}
				}
				iRow[total] = rowindex;
				iCol[total] = colindex;
				//printf("HRow:%d Hcol:%d\n", iRow[total], iCol[total]);
				total++;
			}

		}
	} else {
		//printf("Inside value\n");
		int index = getNumConstraint() - getNumFifthConst();
		int t = 0;
		int id = mtlb->getTargetsId(t);
		int *numvoxel = mtlb->getNumVoxel();
		int numVox = numvoxel[id];
		int totItr = mtlb->getHMatrix()->getNzEntry();

		for (int i = 0; i < mtlb->getHMatrix()->getNzEntry(); i++) {
			//for (int i = 0; i < 280593; i++) {

			if (!step4)
				assert(total <= getNumHess());

			//int rowindex = mtlb->getHMatrix()->rowIndex[i];
			//int colindex = mtlb->getHMatrix()->colIndex[i] - 1;
			int rowindex = (int) mtlb->getHMatrix()->getRowNumber(i);
			int colindex = (int) mtlb->getHMatrix()->getColNumber(i);

			//printf("row:%d col:%d ",rowindex,colindex);
			if (rowindex < NUMBEAMLET) {
				val = 0;
				for (int j = 0; j < mtlb->getNumTargets(); j++) {
					//value[total] = mtlb->getHMatrix()->getValAt(i);
					//printf("**************************************************************************************\n");
					hstime = time(NULL);

					double v = mtlb->getHessian(j)->getAt(rowindex, colindex);
					hendtime = time(NULL);
					double d = difftime(hendtime, hstime);
					totalhash += d;

					//printf("v:%lf  ",v);
					val += lambda[index + j] * v;
					//val += lambda[index+j] * 1;
					//printf("lambda[%d]:%lf\n",(index+j),lambda[index+j]);

				}
				//printf("val:%lf\n",val);
				if (step4 && (rowindex == colindex)) {

					value[rowindex] += val;
					//cout<<"value"<<rowindex<<"]:"<<value[rowindex]<<endl;
				} else {
					value[total] = val;
					total++;
				}
			}

			//printf("i:%d\n",i);

//			if (rowindex >= numVox) {
//				index++;
//				t++;
//				if (t < mtlb->getNumTargets()) {
//					id = mtlb->getTargetsId(t);
//					numVox += numvoxel[id];
//				} else {
////					printf("rowIndex:%d\n",rowindex);
////					printf("i:%d\n",i);
////					printf("nzEntry:%d\n",mtlb->getHMatrix()->getNzEntry());
//					assert(i == mtlb->getHMatrix()->getNzEntry());
//				}
//			}

		}
	}
//	printf("Done Hessian\n");
//printf("total in hessian:%d\n", total);
//printf("Total in hessian:%d\n",total);

	if (!step4)
		assert(total==getNumHess());
//	printf("Elapsed time in Hessian :%lf seconds\n", diff);
//	printf("Elapsed time in Hash :%lf seconds\n", totalhash);
	return true;
}
double StepIImpl::getMialphaMax()
{
	if(miAlpha!=0)
		return miAlpha;
	int *stepIIid = mtlb->getStepIId();
	PrioritizedMain::Vector **D =
				new PrioritizedMain::Vector*[mtlb->getNumStepIIOar()];

	mtlb->getDose(D, w2, stepIIid, mtlb->getNumStepIIOar());
	double **alpha = mtlb->getAlpha();
	int len = D[0]->size() * alpha[0][0];
	long double sum=0;

	D[0]->sort();
	//D[0]->print();

	for(int i=0;i<len;i++)
	{
		int j = D[0]->size()-i-1;
		//printf("va:%lf\n",D[0]->getAt(j));
		sum += D[0]->getAt(j);

	}
	//printf("size:%d\n",D[0]->size());
	//printf("sum:%llf\n",sum);
	//printf("len:%d\n",len);
	miAlpha = sum/len;
	return miAlpha;

}

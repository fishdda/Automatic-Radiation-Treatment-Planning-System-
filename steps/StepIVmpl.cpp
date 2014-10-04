/*
 * StepIVmpl.cpp
 *
 *  Created on: Sep 15, 2012
 *      Author: user
 */

#include "StepIVmpl.h"
#ifdef dummy
StepIVmpl::StepIVmpl(oneVoxel *m, StepIIImpl *stpp3)
#else
		StepIVmpl::StepIVmpl(MatlabCommunicationManager *m, StepIIImpl *stpp3)
#endif
		{
//StepIVmpl::StepIVmpl(oneVoxel *m, StepIImpl *stp2) {
	// TODO Auto-generated constructor stub
	mtlb = m;
	NUMBEAMLET = mtlb->getColInfMat();
	this->step3 = stpp3;
//	printf("StepIVmpl col:%d\n",mtlb->getColInfMat());
//	printf("Numbeamlet:%d\n", NUMBEAMLET);
	printf("Hi from stepIV\n");
	numJacFromNewConst = -1;
	newconsjac=NULL;
}

StepIVmpl::~StepIVmpl() {
	// TODO Auto-generated destructor stub

}
void StepIVmpl::getStartingPoint(double *x,double *z_l,double *z_u) {
	
	//double *ZL = step1->getZL();
	//double *ZU = step1->getZU();
	int nv = getNumVar();
	memcpy(x,step3->xn,sizeof(double)*nv);
	
	//int numconst = getNumConstraint();
	/*for(int i=0;i<numconst;i++)
	{
		//printf("zl[%d]:%lf\n",i,ZL[i]);
		//printf("zu[%d]:%lf\n",i,ZU[i]);
		z_l[i] = 1;
		z_u[i] = 1;
	}*/
}
/**
 * The objective is f(x) = sum{i in RIII} D(w)i.
 *
 */
long double StepIVmpl::evaluateObjective(const double *x) {

	double s = 0;
	//printf("Numbeamlet:%d\n",NUMBEAMLET);
	for (int i = 0; i < NUMBEAMLET; i++) {
		assert(i<getNumVar());
		double m = x[i];
		m = m * m;
		s += m;
	}

	//PrioritizedMain::Vector *w = new PrioritizedMain::Vector(x, NUMBEAMLET);
	//w->print();
	//w = w->square();
	//printf("After\n");
	//w->print();
	//long double s = w->sum();
	//printf("Inside evaluate objective\n");
	//w->print();
	//printf("s:%llf\n",s);
	return s;
	//return step2->evaluateObjective(x);

}
/**
 * Evaluate the gradient of function.
 * The function F(x) = 1/|Vi| *sum{i in RIII} D(w)i
 * The gradient of a function is given by: \/f(x) = (df/dx1, df/dx2,df/dx3...df/dxn )
 * The gradient of the function equals to the sum of the influence matrix column wise  .
 *
 */
void StepIVmpl::evaluategradf(const double *x, double *g) {

//	for(int i=0;i<getNumVar();i++)
//	{
//		g[i]=0;
//	}

	for (int i = 0; i < getNumVar(); i++) {
		if (i < NUMBEAMLET) {
			g[i] = 2 * x[i];
		} else
			g[i] = 0;
	}
	//printf("x[0]:%lf\n",x[0]);
	//PrioritizedMain::Vector *w = new PrioritizedMain::Vector(x, NUMBEAMLET);
	//w->print();
	//w = (*w) * 2;
	//printf("After\n");
	//w->print();
	//printf("Inside gradf\n");
	//w->print();
//	for (int i = 0; i < NUMBEAMLET; i++){
//		//g[i] = w->getAt(i);
//		g[i]=1;
//	}
	//return step2->evaluategradf(x,g);

}

/*
 * Get the number of variables in the optimization.
 */
int StepIVmpl::getNumVar() {

	return step3->getNumVar();

}

/*
 * Get the number of constraints.
 *
 */
int StepIVmpl::getNumConstraint() {

	int m = step3->getNumConstraint();
	//return m;
	return m + this->getNumNewConst();

}

/*
 * @Input: None
 * @Return:
 * 	Number of non-zero entries in Jacobian matrix.
 * @Author:
 * 	Paras Babu Tiwari
 */
int StepIVmpl::getNumJac() {

	int m = step3->getNumJac();
	int n = getNumNewConstJac();
	int total = m+n;
	//printf("Total Jacobian:%d\n",total);
	return total;

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
void StepIVmpl::getBounds(double *x_l, double *x_u, double *g_l, double *g_u) {

	//double moh = step2->getObjectiveFnValue();
	int index = 0;

	//Set the upper and lower bound of infinity or max beamlet weight and 0 for variables.
	double slip = (1 + Constant::SLIP) * (1 + Constant::SLIP)
			* (1 + Constant::SLIP);
	printf("s4:%lf\n",slip);
	step3->getBounds(x_l, x_u, g_l, g_u, slip,Constant::S2);
	index = step3->getNumConstraint();
	PrioritizedMain::Vector **D;
	int num = mtlb->getNumStepIIIOar();
	D = new PrioritizedMain::Vector *[num];
	int *id = mtlb->getStepIIId();
	mtlb->getDose(D, step3->getW3(), id, num);
	//step3->divideByVoxel(D);
	int *numVoxel = mtlb->getNumVoxel();
	for (int i = 0; i < num; i++) {
		int n = numVoxel[id[i]];
		long double s = D[i]->sum();
		long double ub = s/n;

		g_l[index] = -2e19;
		g_u[index] = ub;
		//printf("dose of %d:%lf\n",id[i],g_u[index]);
		index++;
	}

	for (int i = 0; i < num; i++) {
		delete D[i];
	}
	delete[] D;
	assert(index==getNumConstraint());
}
/*
 * Evaluate the constraint.
 * @Param:
 * 	 x: Beamlet weight.
 * 	 g: Output parameter, holds the constraints value at x.
 * 	@author:
 * 		Paras Babu Tiwari
 */
void StepIVmpl::evalg(const double *x, double *g) {
	int total = 0;

	step3->evalg(x, g, total);
	evalNewConst(x, g, total);
	assert(total == this->getNumConstraint());

}
/**
 * yi^(alpha)+1/(1-alpha)*sum{ j in 1 to |Vi|} pji^(alpha)<=Mialpha^{max}.
 * The # of constraints equal to the # of stepII organs.
 *	@return the # of  constraints
 *	@Author: Paras Babu Tiwari
 */
int StepIVmpl::getNumNewConst() {

	int numStepIII = mtlb->getNumStepIIIOar();

	return numStepIII;
}
/**
 * yi^(alpha)+1/(1-alpha)*sum{ j in 1 to |Vi|} pji^(alpha)<=Mialpha^{max}.
 * The # of jacobian entries because of y equals to the # of organs in step II because yi qppears for each i in RII.
 * The # of jacobian entries because of pji equals to the total number of voxels in the step II organs.
 *	@return the # of  Jacobian entries
 *	@Author: Paras Babu Tiwari
 */
int StepIVmpl::getNumNewConstJac() {

	if (numJacFromNewConst >= 0)
		return numJacFromNewConst;
	int numStepIII = mtlb->getNumStepIIIOar();
	int num = mtlb->getNumStepIIIOar();
	double **voxels = mtlb->getVoxels();
	int *numVoxel = mtlb->getNumVoxel();
	Matrices **m = new Matrices*[num];
	int *id = mtlb->getStepIIId();
	int total = 0;
	int count = 0;
	for (int i = 0; i < num; i++) {

		m[i] = mtlb->influence->getSubMatrix(id[i], voxels[id[i]],
				numVoxel[id[i]]);

		int col = m[i]->getCol();
		int row = m[i]->getRow();
		//printf("col:%d\n",col);
		//printf("Row:%d\n",row);
		for (int j = 0; j < col; j++) {
			//printf("j:%d\n",j);
			//printf("Row:%d\n",m[i]->getRow());
			for (int k = 0; k < row; k++) {
				//printf("k:%d\n",k);
				if (m[i]->getAt(k, j) != 0) {
					count++;
					break;
				}
			}
		}

	}
	for (int i = 0; i < num; i++) {
		delete m[i];
	}
	delete[] m;
	numJacFromNewConst = count;
	return numJacFromNewConst;
}
/**
 * Evaluate the constraint yi^(alpha)+1/(1-alpha)*sum{ j in 1 to |Vi|} pji^(alpha)<=Mialpha^{max}.
 * Input:
 *		x: Beamlet weight
 *		total: index in the constraint array.
 *		g: Constraint array(output variable).
 *	@Author: Paras Babu Tiwari
 */
void StepIVmpl::evalNewConst(const double *x, double *g, int& total) {

	int num = mtlb->getNumStepIIIOar();
	PrioritizedMain::Vector **D;
	long double sm = 0;
	D = new PrioritizedMain::Vector *[num];
	int *id = mtlb->getStepIIId();
	mtlb->getDose(D, x, id, num);
	//step3->divideByVoxel(D);
	int bt = total;
	for (int i = 0; i < num; i++) {
		int n = mtlb->getNumVoxels(id[i]);
		long double s = D[i]->sum();
		s = s / n;

		g[total] = s;
		total++;
	}
	for (int i = 0; i < num; i++) {
		delete D[i];
	}
	delete[] D;
//	for(int i=bt;i<total;i++)
//		printf("g[%d]:%lf\n",i,g[i]);
}
/**
 * Evaluate the Jacobian of yi^(alpha)+1/(1-alpha)*sum{ j in 1 to |Vi|} pji^(alpha)<=Mialpha^{max}.
 * The Jacobian of the constraint is same as the gradient of the objective function in the second step.
 * The function F(x) = yi(alpha) + 1/(1-alpha)*|vi| sum{j=1 to |vi|} pji(alpha).
 * The Jacobian of a function is given by: \/f(x) = (df/dx1, df/dx2,df/dx3...df/dxn df/dy1 df/dy2..df/dp1 df/dp2...df/dz1 df/dz2..)
 * The jacobian w.r.to x's are zero as there is no x term in the objective function.
 * The jacobian w.r.to y's equals to one as y's terms are linear.
 * The jacobian w.r.to p's equals to 1/(1-alpha)*|vi| because p's are also linear.
 * The gradient w.r.to z's are zero as there is no z term in the objective function.
 */
void StepIVmpl::handleNewConst(int *iRow, int *iCol, double *value,
		const double *x, int &total, int &l) {
	int num = mtlb->getNumStepIIIOar();

	double **voxels = mtlb->getVoxels();
	int *numVoxel = mtlb->getNumVoxel();
	Matrices **m = new Matrices*[num];
	int *id = mtlb->getStepIIId();
	for (int i = 0; i < num; i++) {

		m[i] = mtlb->influence->getSubMatrix(id[i], voxels[id[i]],
				numVoxel[id[i]]);
		//printf("Influene matrxi of stepIII oar\n");
		//m[i]->printNonZeroEntry();
		//printf("Number of non zero entry of %d:%d\n",id[i],numVoxel[id[i]]);
	}
	//printf("m[i]->getTotalNonZero():%d\n",m[0]->getTotalNonZero());
	bool nonzero = false;
	if (value == NULL) {
		for (int i = 0; i < num; i++) {
			for (int j = 0; j < m[i]->getCol(); j++) {
				for (int k = 0; k < m[i]->getRow(); k++) {
					if (m[i]->getAt(k, j) != 0) {
						iRow[total] = l;
						iCol[total] = j;
						//printf("iRow[%d]:%d %d\n",total,iRow[total],iCol[total]);

						nonzero = true;
						total++;
						break;
					}
				}
				//l++;
			}
			//if(nonzero)
			l++;
			nonzero = false;

		}
	} else {
		int oldtotal = total;
		int nj = getNumNewConstJac();
		if(newconsjac==NULL)
		{			
			newconsjac = new double[nj];
			for (int i = 0; i < num; i++) {
				int n = numVoxel[id[i]];	
				m[i]->colSums(value, total,n);
				
			}
			for(int i=0;i<nj;i++)
			{
				newconsjac[i] = value[oldtotal];
				oldtotal++;
			}
		}
		else
		{
			for(int i=0;i<nj;i++)
			{
				value[total]=newconsjac[i];
				total++;
			}
		}
		
	}
	for (int i = 0; i < num; i++) {

		delete m[i];

	}
	delete[] m;
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
void StepIVmpl::getJacobian(int *iRow, int *iCol, double *value,
		const double *x) {
	int total = 0;
	int l = 0;

	step3->getJacobian(iRow, iCol, value, x, total, l);
	//printf("l:%d\n",l);
	//printf("step3->getNumConstraint():%d\n",step3->getNumConstraint());
	//assert(l == (step3->getNumConstraint()-1));
	this->handleNewConst(iRow, iCol, value, x, total, l);
	//printf("l:%d\n",l);
	assert(l<=this->getNumConstraint());
	//printf("Total:%d\n",total);
	assert(total == this->getNumJac());

}
void StepIVmpl::printDose() {

	//w4->printAll();
	printf("Step1 Dose\n");
	step3->Step2()->Step1()->printDose(w4);
	printf("Step2 Dose\n");
	step3->Step2()->printDose(w4);
	printf("Step3 Dose\n");
	step3->printDose(w4);

	double *metrics = mtlb->getMetrics(step3->Step2()->Step1()->getW1(),
			step3->Step2()->getW2(), step3->getW3(), w4);
	int k = 0;
	for (int i = 0; i < mtlb->getNumTargets(); i++) {
		//printf("D95 of Target %d:%lf\n", i, metrics[i]);
		k++;
	}
	for (int i = 0; i < mtlb->getNumStepIIOar(); i++) {
		//printf("moh of BrainStem:%lf\n", metrics[k]);
		k++;
	}

}
void StepIVmpl::setW(const double *x,const double *z_l,const double *z_u,const double *lambda) {
	w4 = new PrioritizedMain::Vector(x, NUMBEAMLET);
}
/**
 * Get the Hessian. We have stored the hessian in matlab index format, so need to substract one.
 * Hessian of this step is same as that of the step II
 */
bool StepIVmpl::getHessian(int *iRow, int *iCol, double *value,
		double obj_factor, const double *lambda) {

	int total = 0;
	if (value == NULL) {
		for (int i = 0; i < NUMBEAMLET; i++) {

			iRow[total] = i;
			iCol[total] = i;
			//cout<<"irow["<<total<<"]:"<<iRow[total]<<endl;
			//cout<<"icol["<<total<<"]:"<<iCol[total]<<endl;
			total++;

		}
		step3->getHessian(iRow, iCol, value, obj_factor, lambda, total, true);

	} else {

		for (int i = 0; i < NUMBEAMLET; i++) {
			//cout<<value[total]<<endl;
			value[total] = 2 * obj_factor;
			//cout<<"value["<<total<<"]:"<<value[total]<<endl;
			total++;

		}
		step3->getHessian(iRow, iCol, value, obj_factor, lambda, total, true);

	}
	assert(total==getNumHess());
	return true;

}
/*
 * Get number of element in the hessian matrix.
 * @return number of non-zero element in the hessian Matrices
 */
int StepIVmpl::getNumHess() {

	int count = 0;//Number of non-zero element that could be added from step3->getNumHess() function.
	for (int i = 0; i < NUMBEAMLET; i++) {
		double v = mtlb->getHMatrix()->getAt(i, i);
		if (v != 0)
			count++;
	}
	//printf("count:%d\n",count);
	//printf("Numbeamlet:%d\n",NUMBEAMLET);
	int n = (NUMBEAMLET - count) + step3->getNumHess();
	//n = NUMBEAMLET;
	return n;

}

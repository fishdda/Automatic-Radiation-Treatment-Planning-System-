/*
 * OptFormImpl.cpp *
 *  Created on: Nov 20, 2011
 *      Author: user
 */

#include "StepImpl.h"

#ifdef dummy
StepImpl::StepImpl(oneVoxel *m)
#else
StepImpl::StepImpl(MatlabCommunicationManager *m)
#endif
{
	//StepImpl::StepImpl(oneVoxel *m) {
	//StepImpl::StepImpl() :
	//		DummyMatlabComm() {
	//StepImpl::StepImpl() :
	//		oneVoxel() {
	// TODO Auto-generated constructor stub
	//time_t start, end;
	//start = time(NULL);
	//printf("%ld hours since January 1, 1970", start);

	//end = time(NULL);
	//printf("%ld hours since January 1, 1970", end);
	//double diff = difftime(end, start);
	//printf("Time needed for influence Matrices :%lf seconds\n", diff);
	mtlb = m;
	numPBs = mtlb->influence->getCol();
	numTargets = mtlb->getNumTargets();
	targetInf = new Matrices*[numTargets];
	int *numvoxel = mtlb->getNumVoxel();
	for (int i = 0; i < numTargets; i++) {
		int id = mtlb->getTargetsId(i);
		int n = numvoxel[id];

		targetInf[i] = mtlb->influence->getSubMatrix(id, mtlb->voxels[id],
				numvoxel[mtlb->getTargetsId(i)]);
	}
	//Dose = NULL;
	NUMBEAMLET = mtlb->getColInfMat();
	minDose = new double[mtlb->getNumTargets()];
	maxDose = new double[mtlb->getNumTargets()];
	for (int i = 0; i < mtlb->getNumTargets(); i++) {
		minDose[i] = -1;
		maxDose[i] = -1;
	}
	GS = NULL;
	//	printf("**********************************************\n");
	//	int id = mtlb->getTargetsId(0);
	//	int *numvoxel = m->getNumVoxel();
	//	int n = numvoxel[id];
	//	printf("id:%d\n", id);
	//	printf("n:%d\n", n);
	//	double *row = m->getVoxels()[id];
	//	for (int t = 0; t < n; t++) {
	//		printf("%lf\n", row[t]);
	//	}
}

StepImpl::~StepImpl() {
	// TODO Auto-generated destructor stub
}

long double StepImpl::evaluateObjective(const double *x) {
	long double sm = 0;

	//printf("Inside evaluate objective\n");
	//time_t start_t, end_t;
	//start_t = time(NULL);
	//      FileManager *f = new FileManager("/Users/user/CoinIpopt/Ipopt/examples/Project/Debug/w.txt");
	////    double *result = new double[influence->getCol()];
	////    f->readFile(result);
	//Matrices *xm = new Matrices(x, mtlb->influence->getCol());
	int numTarget = mtlb->getNumTargets();
	PrioritizedMain::Vector **D = new PrioritizedMain::Vector*[numTarget];
	getTargetDose(x, D);
	this->subPrescDose(D);
	//printf("Before square\n");
	//D[0]->print();
	for (int i = 0; i < numTarget; i++) {
		D[i] = D[i]->square();
		//D[i] * D[i];
	}
	//printf("After square\n");
	//D[0]->print();
	//this->divideByVoxel(D);
	int *numvoxel = mtlb->getNumVoxel();
	double *prescDose = mtlb->getPrescribedDose();
	int numPBs = mtlb->influence->getCol();

	for (int i = 0; i < numTarget; i++) {
		// D[i]->print();
		long double s = D[i]->sum();
		int n = numvoxel[mtlb->getTargetsId(i)];
		s = s / n;
		sm += s;
		//printf("******* Sum:%lf ***********\n",sm);
	}

	double t = 0;
	for (int i = 0; i < mtlb->getNumTargets(); i++) {
		//t += 1.0 / 4 * x[numPBs + i] * x[numPBs + i];
		t += x[numPBs + i] * x[numPBs + i];
		//t += x[numPBs + i] ;
	}
	sm = sm + t;
	for (int i = 0; i < numTarget; i++) {
		delete D[i];
	}
	delete[] D;
	
	//delete D;
	//delete xm;
	//end_t = time(NULL);
	//double diff = difftime(end_t, start_t);
	//printf("objective elapsed time:%lf\n", diff);
	return sm;

}
/*
 * This function evaluates Gi(w) = 1/|Vi| sum{j in Vi}[Dj(w)-Dipre]^2. The w matrix is obtained through optimization.
 * @param: Output Array G.
 * @return Gi(w) array.
 *
 */
void StepImpl::evaluateG(double *G, PrioritizedMain::Vector *w) {
	int numTarget = mtlb->getNumTargets();
	PrioritizedMain::Vector **D = new PrioritizedMain::Vector*[numTarget];
	//Matrices *D = new Matrices(mtlb->getTTLTargetVoxel(), 1);
	assert(G != NULL);
	getTargetDose(w, D);

	//printf("Max Dose:%lf",D[0]->getMax());
	//printf("Average Dose:%lf",D[0]->average());
	//printf("Min Dose:%lf",D[0]->getMin());
	int start = 0, end = 0, numVoxel = 0;
	this->subPrescDose(D);

	for (int i = 0; i < numTarget; i++) {

		D[i] = D[i]->square();	//->vectorMultiply(D, D, D);

	}
//	printf("D after multiply\n");
//	D->printNonZeroEntry();
	this->divideByVoxel(D);

//	printf("D after divide\n");
//	D->printNonZeroEntry();
	for (int i = 0; i < numTarget; i++) {
		G[i] = D[i]->sum();

	}
	for (int i = 0; i < numTarget; i++) {
		delete D[i];
	}
	delete[] D;
	//this->sumByTarget(D, G);
	//delete D;
}
/*
 * This function evaluates Gi(w) = 1/|Vi| sum{j in Vi}[Dj(w)-Dipre]^2. The w matrix is obtained through optimization.
 * @param: Output Array G.
 * @return Gi(w) array.
 *
 */
 /*
void StepImpl::evaluateG(double*& G, PrioritizedMain::Vector *w) {

	if (GS != NULL) {
		//printf("copying\n");
		G = GS;
	} else {
		int numTarget = mtlb->getNumTargets();
		PrioritizedMain::Vector **D = new PrioritizedMain::Vector*[numTarget];
		//Matrices *D = new Matrices(mtlb->getTTLTargetVoxel(), 1);
		//assert(G != NULL);
		//printf("c1\n");
		GS = new double[numTarget];
		getTargetDose(w, D);
		//printf("c2\n");
		//printf("Max Dose:%lf",D[0]->getMax());
		//printf("Average Dose:%lf",D[0]->average());
		//printf("Min Dose:%lf",D[0]->getMin());
		int start = 0, end = 0, numVoxel = 0;
		this->subPrescDose(D);

		for (int i = 0; i < numTarget; i++) {

			D[i] = D[i]->square(); //->vectorMultiply(D, D, D);

		}
		//	printf("D after multiply\n");
		//	D->printNonZeroEntry();
		this->divideByVoxel(D);
		//printf("c2.1\n");
		//	printf("D after divide\n");
		//	D->printNonZeroEntry();
		for (int i = 0; i < numTarget; i++) {
			GS[i] = D[i]->sum();

		}
		//printf("c3\n");
		for (int i = 0; i < numTarget; i++) {
			delete D[i];
		}
		delete[] D;
		G = GS;
	}

	//this->sumByTarget(D, G);
	//delete D;
	//printf("Done evalG\n");
}*/
/**
 * Sum the dose of a target
 * @param: D Matrices: Dose of individual voxel. The matrix must be of Total # of voxels in all targets * 1
 * @param: Output Array G.
 * @return: An array of size equals to the # of targets containing dose of each target.
 */
//void StepImpl::sumByTarget(Matrices *D, double *G) {
//	int start = 0, end = 0, numVoxel = 0;
//	int *numvoxel = mtlb->getNumVoxel();
//	for (int i = 0; i < mtlb->getNumTargets(); i++) {
//		start += numVoxel;
//		numVoxel = numvoxel[mtlb->getTargetsId(i)];
//		end = start + numVoxel;
//		G[i] = D->sum(start, end);
//
//	}
//}
/*
 * Divide the matrix D by number of voxel.
 * The D corresponds to the dose to different target structures.
 * So different target structure might have different number of voxel.
 * The code finds the number of voxel corresponding to dose entry and divide it by that number.
 */
void StepImpl::divideByVoxel(PrioritizedMain::Vector **D) {
	int start = 0, end = 0, n = 0;
	int *numvoxel = mtlb->getNumVoxel();
	for (int i = 0; i < mtlb->getNumTargets(); i++) {
		//start += numVoxel;
		n = numvoxel[mtlb->getTargetsId(i)];
		//end = start + numVoxel;
		D[i] = (*D[i]) / n; //->divide(numVoxel, start, end);
		//printf("n:%d\n",n);

	}
}
/*
 * Get number of element in the hessian matrix.
 * @return number of non-zero element in the hessian Matrices
 */
int StepImpl::getNumHess() {
	//return 0;//limited memory;

	if (mtlb->getHMatrix() == NULL)
		return 0;
	//For real data
	//
#ifdef dummy
	return mtlb->getHMatrix()->getNzEntry() + mtlb->getNumTargets();
#else
	return mtlb->getHMatrix()->getNzEntry();
#endif

}

/**
 * Get the Hessian. We have stored the hessian in matlab index format, so need to substract one.
 * Hessian is A'*A
 */
bool StepImpl::getHessian(int *iRow, int *iCol, double *value,
		double obj_factor, const double *lambda) {

	//return false;//limited memory
	//time_t start_t, end_t;
	//start_t = time(NULL);
	int total = 0;
	if (mtlb->getHMatrix() == NULL)
		return false;
	if (value == NULL) {
		for (int i = 0; i < mtlb->getHMatrix()->getNzEntry(); i++) {
			assert(total < getNumHess());
			//			iRow[total] = mtlb->getHMatrix()->rowIndex[i] - 1;
			//			iCol[total] = mtlb->getHMatrix()->colIndex[i] - 1;
			iRow[total] = (int) mtlb->getHMatrix()->getRowNumber(i);
			iCol[total] = (int) mtlb->getHMatrix()->getColNumber(i);
			//printf("row:%d col:%d\n", iRow[total], iCol[total]);
			total++;
		}
#ifdef dummy
		//Comment out this when using the real data
		for (int j = 0; j < mtlb->getNumTargets(); j++) {
			iRow[total] = mtlb->influence->getCol() + j;
			iCol[total] = mtlb->influence->getCol() + j;
			total++;
		}
#endif
	} else {
		for (int i = 0; i < mtlb->getHMatrix()->getNzEntry(); i++) {
			assert(total < getNumHess());
			value[total] = obj_factor * mtlb->getHMatrix()->getValAt(i);
			//printf("value[%d]:%lf\n",total,value[total]);
			total++;

		}
#ifdef dummy

		//Comment out this when using the real data
		for (int j = 0; j < mtlb->getNumTargets(); j++) {
			value[total] = obj_factor * 0.5;

			total++;
		}
#endif
	}

	//end_t = time(NULL);
	//double diff = difftime(end_t, start_t);
	//printf("Time in Hessian:%lf\n", diff);
	//printf("Total in hessian:%d\n",total);
	assert(total == getNumHess());
	return true;
}
void StepImpl::getStartingPoint(double *x,double *z_l,double *z_u) {
//void StepImpl::getStartingPoint(double *x) {
	for (int i = 0; i < getNumVar(); i++) {
		x[i] = 2;
		
	}
}
void StepImpl::setW(const double *x,const double *z_l,const double *z_u,const double *lambda)
{
	w1 = new PrioritizedMain::Vector(x, NUMBEAMLET);
	int numconst = getNumConstraint();
	/*
	ZL = new double[numconst];
	ZU = new double[numconst];
	for(int i=0;i<numconst;i++)
	{
		printf("zl[%d]:%lf\n",i,z_l[i]);
		printf("zu[%d]:%lf\n",i,z_u[i]);
		ZL[i] = z_l[i];
		ZU[i] = z_u[i];
	}*/
	
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
void StepImpl::getBounds(double *x_l, double *x_u, double *g_l, double *g_u) {
	for (int i = 0; i < mtlb->getColInfMat(); i++) {
		x_l[i] = 0.0;
		x_u[i] = mtlb->getUBWeight();
	}
	double *prescDose = mtlb->getPrescribedDose();
	for (int i = 0; i < mtlb->getNumTargets(); i++) {
		int ix = mtlb->getColInfMat() + i;
		x_l[ix] = 0.05 * prescDose[i];
		//x_l[ix] = 0;
		//x_u[ix] = mtlb->getUBWeight();
		x_u[ix] = 2e19;
	}
	// the first constraint g1 has a lower bound of infinity
	// The upper bound for the first constraint is 0
	for (int i = 0; i < mtlb->getTTLTargetVoxel(); i++) {
		g_l[i] = -2e19;
		g_u[i] = 0;

	}

	int no = mtlb->getTTLTargetVoxel();
	// the second constraint g2 is an equality constraint, so we set the
	// upper and lower bound to the same value
	int maxItr = no + mtlb->getTTLStepIOar();
	for (int i = no; i < maxItr; i++) {
		g_l[i] = -2e19;
		g_u[i] = mtlb->getMaxDose();
		//printf("Max Dose:%lf\n",g_u[i]);
	}
}
/*
void StepImpl::getJacobian(int *iRow, int *iCol, double *value, const double *x) {

	int numPBs = mtlb->influence->getCol();
	int total = 0;
	//time_t start_t, end_t;
	//start_t = time(NULL);
	double val = 0;
	float v = 0;
	int l = 0;
	
	for (int i = 0; i < numTargets; i++) 
	{
	
		if(value==NULL)
		{
				std::vector<long int> rows = mtlb->targetInf[i]->getNzRows();
				std::vector<long int> cols = mtlb->targetInf[i]->getNzCols();
				int sz = mtlb->targetInf[i]->value.size();
				int numrow = mtlb->targetInf[i]->row;
				int prevrow = 0;
				for (int k = 0; k < sz; k++) 
				{
					if(prevrow!=0)
					{
						if(prevRow != rows[k])
						{
							iRow[total] = l + prevRow;
							iCol[total] = numPBs+1;
							mtlb->targetInf[i]->negValue.push_back(-1);
							total++;
						}
					}
					iRow[total] = l + rows[k];
					iCol[total] = cols[k];					
					prevrow = rows[k];
					negValue.push_back(mtlb->targetInf[i]->value[k]);
					total++;
				}
				l = l + numrow;
				
		}
		else
		{
			int sz = mtlb->targetInf[i]->value.size() + numrow;
			memcpy(value + total, negValue, sizeof(double) * sz);
			total = total + sz;
			l = l + mtlb->targetInf[i]->row;
		}		
	}
	int numStepIOar = mtlb->getNumStepIOar();
	for (int i = 0; i < numStepIOar; i++) 
	{
		
		if(value==NULL)
		{
			std::vector<long int> rows = mtlb->stepIinf[i]->getNzRows();
			std::vector<long int> cols = mtlb->stepIinf[i]->getNzCols();
			int sz = mtlb->stepIinf[i]->value.size();
			int numrow = mtlb->stepIinf[i]->row;
			for (int k = 0; k < sz; k++) 
			{
				iRow[total] = l + rows[k];
				iCol[total] = cols[k];	
				mtlb->stepIinf[i]->negValue.push_back(mtlb->stepIinf[i]->value[k]);				
				total++;
			}
			l = l + numrow;
		}
		else
		{
			int sz = mtlb->targetInf[i]->value.size() ;
			memcpy(value + total, negValue, sizeof(double) * sz);
			total = total + sz;
			l = l + mtlb->stepIinf[i]->row;
		}
	}
	assert(total == getNumJac());

}*/
/**
 *
 */
void StepImpl::getJacobian(int *iRow, int *iCol, double *value, const double *x) {

	int numPBs = mtlb->influence->getCol();
	int total = 0;
	//time_t start_t, end_t;
	//start_t = time(NULL);
	double val = 0;
	float v = 0;
	int l = 0;
	if (value == NULL) {
		int *numvoxel = mtlb->getNumVoxel();
		for (int i = 0; i < mtlb->getNumTargets(); i++) {

			int id = mtlb->getTargetsId(i);
			int n = numvoxel[id];
			//			printf("idL%d\n", id);
			//			printf("n:%d\n", n);
			double *row = mtlb->getVoxels()[id];
			//printf("id:%d\n",id);
			//printf("Row address:%x\n",row);
			//			for (int t = 0; t < n; t++) {
			//				printf("%lf\n", row[t]);
			//			}
			//printf("Number of Voxel %d of target:%d\n", n, id);

			for (int k = 0; k < n; k++) {
				//printf("r before:%lf\n",row[k]);
				long unsigned int r = (long unsigned int) row[k];
				//printf("r:%lu\n", r);
				for (int j = 0; j < numPBs; j++) {

					assert(r >= 0);
					val = mtlb->influence->getAt(r, j);

					if (val != 0) {
						///printf("r:%d col:%d:%lf\n",r,j,val);
						//printf("Row:%d,col:%d:%lf\n",r,j,val);
						iRow[total] = l; //r;
						iCol[total] = j; //j;
						total++;
						//						printf("Total:%d\n",total);
						//						printf("row:%d col:%d:%lf\n",r,j,val);

					}

				}
				iRow[total] = l;
				iCol[total] = numPBs + i;
				total++;
				//printf("Total:%d\n",total);
				l++;

			}
		}
		//	printf("total:%d\n",total);
		l = 0;
		//for (int i = 0; i < mtlb->getTTLStepIOar(); i++) {
		for (int i = 0; i < mtlb->getNumStepIOar(); i++) {
			//int r = getRowStepI(i);
			//assert(r >= 0);
			int id = mtlb->getStepIOrganId(i);
			double *row = mtlb->getVoxels()[id];
			int n = numvoxel[id];
			for (int k = 0; k < n; k++) {
				long unsigned int r = row[k];
				for (int j = 0; j < numPBs; j++) {

					val = mtlb->influence->getAt(r, j);

					if (val != 0) {
						//printf(" stepI row:%d col:%d:%lf\n",r,j,val);
						iRow[total] = l + mtlb->getTTLTargetVoxel(); //r;
						iCol[total] = j;
						total++;

					}

				}
				l++;
			}
		}

	} else {
		int *numvoxel = mtlb->getNumVoxel();
		for (int i = 0; i < mtlb->getNumTargets(); i++) {
			int id = mtlb->getTargetsId(i);
			int n = numvoxel[id];
			double *row = mtlb->getVoxels()[id];
			for (int k = 0; k < n; k++) {
				long unsigned int r = row[k];
				for (int j = 0; j < numPBs; j++) {
					val = mtlb->influence->getAt(r, j);
					if (val != 0) {
						value[total] = -val;
						total++;
					}

				}
				value[total] = -1;
				total++;

			}
		}
		for (int i = 0; i < mtlb->getNumStepIOar(); i++) {
			int id = mtlb->getStepIOrganId(i);
			double *row = mtlb->getVoxels()[id];
			int n = numvoxel[id];

			for (int k = 0; k < n; k++) {
				//for (int i = 0; i < mtlb->getTTLStepIOar(); i++) {
				//int r = getRowStepI(i);
				long unsigned int r = row[k];
				for (int j = 0; j < numPBs; j++) {
					val = mtlb->influence->getAt(r, j);

					if (val != 0) {
						v = val;
						//printf("v:%f,val:%lf\n",v,val);
						value[total] = val;
						//value[total] = v;

						total++;
					}

				}
			}
		}
	}
	//end_t = time(NULL);
	//double diff = difftime(end_t, start_t);
	//printf("Time in Jacobian:%lf\n", diff);
	//printf("total:%d\n", total);
	//printf("getNumJac():%d\n",getNumJac());
	assert(total == getNumJac());

} 

void StepImpl::evalg(const double *x, double *g) {
	int s;

	assert(g != NULL);
	assert(x != NULL);
	time_t start_t, end_t;
	//start_t = time(NULL);
	int numTarget = mtlb->getNumTargets();
	s = mtlb->getColInfMat() + numTarget;

	int numPBs = mtlb->influence->getCol();
	PrioritizedMain::Vector **D = new PrioritizedMain::Vector *[numTarget];
	this->getTargetDose(x, D);
	//end_t = time(NULL);
	double diff = difftime(end_t, start_t);
	//printf("target dose elapsed time:%lf\n", diff);
	
	//start_t = time(NULL);
	int numVoxel = 0;
	//this->subPrescDose(D);
	double *prescDose = mtlb->getPrescribedDose();
	for (int i = 0; i < mtlb->getNumTargets(); i++) {
		D[i] = D[i]->subFrom(prescDose[i]);
	}
	//end_t = time(NULL);
	//diff = difftime(end_t, start_t);
	//printf("sub dose elapsed time:%lf\n", diff);
	
	
	int *numvoxel = mtlb->getNumVoxel();
	int l = 0;
	//printf("t:%lf\n",x[numPBs+0]);
	for (int i = 0; i < mtlb->getNumTargets(); i++) {
		numVoxel = numvoxel[mtlb->getTargetsId(i)];

		//printf("Before substraction\n");
		//D[i]->print();
		D[i] = (*D[i]) - x[numPBs + i];
		//printf("After substraction\n");
		//D[i]->print();
		for (int k = 0; k < D[i]->size(); k++) {
			g[l] = D[i]->getAt(k);
			l++;
		}
		//printf("After subtracting\n");
		//D[i]->print();

	}
	
	//start_t = time(NULL);
	PrioritizedMain::Vector **DS =
			new PrioritizedMain::Vector *[mtlb->getTTLStepIOar()];
	mtlb->getDose(DS, x, mtlb->getStepIOrganId(), mtlb->getNumStepIOar());
	//int k = 0;
	//end_t = time(NULL);
	//diff = difftime(end_t, start_t);
	//printf("first step dose elapsed time:%lf\n", diff);
	
	//start_t = time(NULL);
	
	int numStepIOar = mtlb->getNumStepIOar();
	for (int i = 0; i < numStepIOar; i++) {
		double *d = DS[i]->getDblArray();
		int s = DS[i]->size();
		memcpy(g+l,d,sizeof(double) * s);
		l = l + s;
		/*for (int k = 0; k < DS[i]->size(); k++) {
			g[l] = DS[i]->getAt(k);
			l++;
		}*/
	}
	//printf("deleting D\n");
	for (int i = 0; i < numTarget; i++) {
		delete D[i];
	}
	delete[] D;
	////	printf("Deleting DS\n");
	//	for (int i = 0; i < mtlb->getNumStepIOar(); i++) {
	//		delete DS[i];
	//	}
	//	delete[] DS;
	assert(l == getNumConstraint());
	//	for (int i = 0; i < numTarget; i++) {
	//		memcpy(&g[k], D[i]->getDblArray(), sizeof(double) * D[i]->size());
	//		k = k + D[i]->size();
	//	}
	//	int numStepIOar = mtlb->getNumStepIOar();
	//	for (int i = 0; i < numStepIOar; i++) {
	//		memcpy(&g[k], DS[i]->getDblArray(), sizeof(double) * DS[i]->size());
	//		k = k + DS[i]->size();
	//	}
		//for (int i = 0; i < getNumConstraint(); i++)
	//		printf("g[%d]:%lf\n", i, g[i]);
	//cout<<"g["<<i<<"]:"<<g[i]<<endl;

	//delete D;
	//delete DS;
	//printf("Done with eval g\n");
	//delete constraint;
	//end_t = time(NULL);
	//diff = difftime(end_t, start_t);
	//printf("copy elapsed time:%lf\n", diff);


}
//Matrices *StepImpl::getwMatrix(const double *x) {
//	return new Matrices(x, mtlb->influence->getCol() + mtlb->getNumTargets());
//}
/*
 * Get the number of voxels in a target.
 * ith voxels of the target. This is not a row number in the influence matrix rather it's the row number in the matrix D.
 * Number of voxels in the target
 */
//int StepImpl::getNumVoxels(int i) {
//	int num = 0;
//	int dividend = 0;
//	int remainder = i;
//	int rowInd = 0;
//	int *numvoxel = mtlb->getNumVoxel();
//	for (int j = 0; j < mtlb->getNumTargets(); j++) {
//		dividend += numvoxel[mtlb->getTargetsId(j)];
//		int result = i / dividend;
//		if (result == 0) {
//			int id = mtlb->getTargetsId(j);
//			num = numvoxel[id];
//			break;
//		}
//	}
//	return num;
//}
/**
 * Get the row number in the influence matrix corresponding to the row number of the D matrix.
 * @param ith row of the D matrix
 * @return row number of the influence matrix
 */
//int StepImpl::getRow(int i) {
//	int dividend = 0;
//	int remainder = i;
//	int rowInd = 0;
//	int *numvoxel = mtlb->getNumVoxel();
//	for (int j = 0; j < mtlb->getNumTargets(); j++) {
//		dividend += numvoxel[mtlb->getTargetsId(j)];
//		int result = i / dividend;
//		if (result == 0) {
//			dividend -= numvoxel[mtlb->getTargetsId(j)];
//			remainder = i - dividend;
//			int id = mtlb->getTargetsId(j);
//			double *row = mtlb->voxels[id];
//			rowInd = (int) row[remainder];
//			break;
//		}
//	}
//
//	return rowInd;
//}
//int StepImpl::getRowStepI(int i) {
//	int dividend = 0;
//	int remainder = i;
//	int rowInd = 0;
//	int *numvoxel = mtlb->getNumVoxel();
//	int *stepIOarId = mtlb->getStepIOrganId();
//	for (int j = 0; j < mtlb->getNumStepIOar(); j++) {
//		int k = (int) stepIOarId[j];
//		dividend += numvoxel[k];
//		int result = i / dividend;
//		if (result == 0) {
//			dividend -= numvoxel[k];
//			remainder = i - dividend;
//
//			double *row = mtlb->voxels[k];
//
//			rowInd = (int) row[remainder];
//
//		}
//	}
//
//	return rowInd;
//}
int StepImpl::getNumJac() {

	//	printf("Total Target Voxel:%d\n",mtlb->getTTLTargetVoxel());
	//	printf("Non zero voxel of target:%d\n", mtlb->getNonZeroVoxelofTargets());
	//	printf("Non zero voxel of stepI:%d\n",mtlb->getNonzeroStepIOarVoxels());
	int nnz_jac_g = mtlb->getTTLTargetVoxel()
			+ mtlb->getNonZeroVoxelofTargets()
			+ mtlb->getNonzeroStepIOarVoxels();
	;
	//int nnz_jac_g = mtlb->getNonZeroVoxelofTargets();
	//+ mtlb->getNonzeroStepIOarVoxels();
	return nnz_jac_g;

}

int StepImpl::getNumVar() {
	//printf("mtlb:%x\n",mtlb);
	//printf("col:%d\n",mtlb->getColInfMat());
	//printf("# of Targets:%d\n",mtlb->getNumTargets());

	return mtlb->getColInfMat() + mtlb->getNumTargets();
	//return this->getColInfMat() ;
}
int StepImpl::getNumConstraint() {
	int m = mtlb->getTTLTargetVoxel() + mtlb->getTTLStepIOar();
	//int m = this->getTTLTargetRow();
	return m;

}
/**
 * Evaluate the gradient of function.
 * Let the influence matrix size is m(total number of voxels in all targets).
 * 1<=i<= numPBs
 * g[i] = 2* D1 * a1i + 2* D2 * a2i + 2* D3 *a3i+...+ Dm*ami
 * It's caller responsibility to allocate the space for g.
 * @Author: Paras Babu Tiwari
 */
void StepImpl::evaluategradf(const double *x, double *g) {

	assert(g != NULL);
	time_t start_t, end_t;
	Matrices *res;
	start_t = time(NULL);
	int numPBs = mtlb->getColInfMat();
	int numTarget = mtlb->getNumTargets();
	PrioritizedMain::Vector **D = new PrioritizedMain::Vector *[numTarget]; //Matrices(mtlb->getTTLTargetVoxel(), 1);
	PrioritizedMain::Vector **grad = new PrioritizedMain::Vector *[numTarget];
	getTargetDose(x, D);
	this->subPrescDose(D);
	this->divideByVoxel(D);
	double coeff = 0;

	int l = 0;
	int n = 0;
	double val = 0, tempVal = 0;
	float v = 0;
	for (int j = 0; j < getNumVar(); j++) {
		assert(j < getNumVar());
		//printf("x[%d]:%lf\n",j,x[j]);
		g[j] = 0;

	}

	int *numvoxel = mtlb->getNumVoxel();
	//start_t = time(NULL);
	for (int i = 0; i < mtlb->getNumTargets(); i++) {
		int id = mtlb->getTargetsId(i);
		res = mtlb->influence->getSubMatrix(id, mtlb->voxels[id],
				numvoxel[mtlb->getTargetsId(i)]);
		//printf("Row:%d Col:%d\n",res->getRow(),res->getCol());

		grad[i] = res->matrixTransMultiply(D[i]);
		//D[i]->print();
		//printf("grad[i]: %d\n",grad[i]->size());
		//res->printNonZeroEntry();
		//		for(int j=0;j<grad[i]->size();j++)
		//		{
		//			printf("grad[%d]:%lf\n",j,grad[i]->getAt(j));
		//		}
		delete res;
	}
	//end_t = time(NULL);
	//double diff = difftime(end_t, start_t);
	//printf("gradient computation time:%lf\n", diff);
	start_t = time(NULL);
	for (int i = 0; i < mtlb->getNumTargets(); i++) {
		for (int j = 0; j < mtlb->influence->getCol(); j++) {
			double v = grad[i]->getAt(j);
			//printf("v:%lf\n",v);
			v = v * 2;

			g[j] = g[j] + v;

		}

	}
	end_t = time(NULL);
	//diff = difftime(end_t, start_t);
	//printf("gradient copy time:%lf\n", diff);
	for (int i = 0; i < mtlb->getNumTargets(); i++) {
		assert((numPBs + i) < getNumVar());
		//g[numPBs + i] = 1.0 / 2 * x[numPBs + i];
		g[numPBs + i] = 2 * x[numPBs + i];
		//g[numPBs + i] =  1;

	}
	for (int i = 0; i < numTarget; i++) {
		delete D[i];
	}
	delete[] D;
	end_t = time(NULL);
	//diff = difftime(end_t, start_t);
	//printf("gradient time:%lf\n", diff);
}

void StepImpl::evaluategradf1(const double *x, double *g) {

	//evaluategradf(x,g);
	assert(g != NULL);
	time_t start_t, end_t;
	start_t = time(NULL);
	int numPBs = mtlb->getColInfMat();
	//Matrices *xm = new Matrices(x, mtlb->influence->getCol());
	int numTarget = mtlb->getNumTargets();
	PrioritizedMain::Vector **D = new PrioritizedMain::Vector *[numTarget]; //Matrices(mtlb->getTTLTargetVoxel(), 1);
	getTargetDose(x, D);
	this->subPrescDose(D);
	this->divideByVoxel(D);
	double coeff = 0;

	int l = 0;
	int n = 0;
	double val = 0, tempVal = 0;
	float v = 0;
	for (int j = 0; j < getNumVar(); j++) {
		assert(j < getNumVar());
		g[j] = 0;

	}
	end_t = time(NULL);
	double diff = difftime(end_t, start_t);
	printf("first time:%lf\n", diff);
	start_t = time(NULL);
	int *numvoxel = mtlb->getNumVoxel();
	for (int i = 0; i < mtlb->getNumTargets(); i++) {
		int id = mtlb->getTargetsId(i);
		n = numvoxel[id];
		double *row = mtlb->getVoxels()[id];
		l = 0;
		for (int k = 0; k < n; k++) {
			int r = row[k];
			for (int j = 0; j < mtlb->influence->getCol(); j++) {
				//for (int j = 0; j < 1100; j++) {
				//		int r = row[k];
				//cout<<"j:"<<j<<endl;
				coeff = mtlb->influence->getAt(r, j);

				if (coeff != 0) {
					//printf("coeff[%d][%d]:%lf\n",r,j,coeff);
					assert(j < getNumVar());
					coeff = coeff * 2;

					val = D[i]->getAt(l) * coeff;
					//val = coeff;


					assert(j < numPBs);
					g[j] += val;

				}
			}

			l++;

		}
		//		for (int j = 0; j < getNumVar(); j++) {
		//			assert(j < getNumVar());
		//			g[j] = g[j] / n;
		//
		//		}


	}
	end_t = time(NULL);
	diff = difftime(end_t, start_t);
	printf("second time:%lf\n", diff);
	for (int i = 0; i < mtlb->getNumTargets(); i++) {
		assert((numPBs + i) < getNumVar());
		//g[numPBs + i] = 1.0 / 2 * x[numPBs + i];
		g[numPBs + i] = 2 * x[numPBs + i];
		//g[numPBs + i] =  1;

	}
	for (int i = 0; i < numTarget; i++) {
		delete D[i];
	}
	delete[] D;
	//delete xm;
	//delete D;
	//end_t = time(NULL);
	//double diff = difftime(end_t, start_t);
	//printf("gradient elapsed time:%lf\n", diff);
}
void StepImpl::print() {
	//	Matrices *m = new Matrices(getNumVoxels(getTargetsId(0)),
	//			influence->getCol());
	//
	//	this->influence->getSubMatrix(voxels[0], getNumVoxels(getTargetsId(0)),
	//			m);
	//	m->printMatrix();
	//	delete m;
}
/*
 * Get the dose of the step I.
 * @param:
 * 	w: Beamlet weight.
 * 	D: Output variable.The dose matrix.
 * @return:None
 */
//void StepImpl::getStepIDose(Matrices *w, Matrices *D) {
//
//	assert(D != NULL);
//	assert(D->getRow() > 0);
//	int offset = 0;
//	Matrices *m;
//	time_t start_t, end_t;
//	start_t = time(NULL);
//	int *numvoxel = mtlb->getNumVoxel();
//	for (int i = 0; i < mtlb->getNumStepIOar(); i++) {
//
//		int id = mtlb->getStepIOrganId(i);
//		m = mtlb->influence->getSubMatrix(id, mtlb->voxels[id],
//				numvoxel[mtlb->getStepIOrganId(i)]);
//		Matrices *result = new Matrices(m->getRow(), w->getCol());
//		result->matrixMultiply(m, w);
//		assert(result->getRow() <= D->getRow());
//		D->copy(result, offset);
//		offset = numvoxel[mtlb->getStepIOrganId(i)];
//		delete result;
//	}
//	end_t = time(NULL);
//	double diff = difftime(end_t, start_t);
//
//}
void StepImpl::getTargetDose(PrioritizedMain::Vector *w,
		PrioritizedMain::Vector **D) {
	assert(D != NULL);
	//assert(D->getRow() > 0);
	struct timeval start, end;
	long mtime, seconds, useconds;
	//int offset = 0;
	//Matrices *res;
	//time_t start_t, end_t;
	//start_t = time(NULL);
	int *numvoxel = mtlb->getNumVoxel();
	//cout<<"w"<<endl;
	//w->print();
	for (int i = 0; i < mtlb->getNumTargets(); i++) {
		//		int id = mtlb->getTargetsId(i);
		//		//	printf("Target id:%d\n",id);
		//		res = mtlb->influence->getSubMatrix(id, mtlb->voxels[id],
		//				numvoxel[mtlb->getTargetsId(i)]);
		//		printf("Res\n");
		//		res->printNonZeroEntry();
		//res->print();
		//assert(result->getRow() <= D->getRow());
		//D->copy(result, offset);
		
		D[i] = targetInf[i]->matrixMultiply(w);
		
		//D[i]->print();
		//offset += numvoxel[mtlb->getTargetsId(i)];
		//delete result;
		//end_t = time(NULL);
		//double diff = difftime(end_t, start_t);
		//delete res;

	}
	//printf("Returned from get dose\n");
}
void StepImpl::getTargetDose(const double *x, PrioritizedMain::Vector **D) {

	assert(D != NULL);
	//assert(D->getRow() > 0);
	//double *xt = new double[mtlb->getColInfMat()];
	//	for(int i=0;i<mtlb->getColInfMat();i++)
	//		xt[i] = i+1;
	PrioritizedMain::Vector *w = new PrioritizedMain::Vector(x,
			mtlb->getColInfMat());
	//w->print();
	//std::cout<<"w size"<<w->V().size()<<endl;
	getTargetDose(w, D);
	delete w;
	//	int offset = 0;
	//	Matrices *res;
	//	time_t start_t, end_t;
	//	start_t = time(NULL);
	//	int *numvoxel = mtlb->getNumVoxel();
	//	for (int i = 0; i < mtlb->getNumTargets(); i++) {
	//		int id = mtlb->getTargetsId(i);
	//		res = mtlb->influence->getSubMatrix(id, mtlb->voxels[id],
	//				numvoxel[mtlb->getTargetsId(i)]);
	////		printf("Res\n");
	////		res->printNonZeroEntry();
	//		PrioritizedMain::Vector *result = new PrioritizedMain::Vector(
	//				res->getRow());
	//		res->matrixMultiply(w, result);
	//		//assert(result->getRow() <= D->getRow());
	//		//D->copy(result, offset);
	//		D[i] = result;
	//		//offset += numvoxel[mtlb->getTargetsId(i)];
	//		//delete result;
	//		//end_t = time(NULL);
	//		//double diff = difftime(end_t, start_t);
	//
	//	}
	//
	//	return D;
}
void StepImpl::subPrescDose(PrioritizedMain::Vector **D) {

	int *numVoxel = mtlb->getNumVoxel();
	double *prescDose = mtlb->getPrescribedDose();
	int start = 0;
	int end = 0;
	for (int i = 0; i < mtlb->getNumTargets(); i++) {
		//		int id = mtlb->getTargetsId(i);
		//		end += numVoxel[id];
		//D[i]->substract(prescDose[i]);
		D[i] = (*D[i]) - prescDose[i];
		//printf("prescribtion dose:%lf\n",prescDose[i]);
		//§D[i] = (*D[i]) * -1;
		//start +=end;

	}
}
/*
 * Print the dose
 */
void StepImpl::printDose() {

	printDose(w1);
	//w1->printAll();
}
/*
 * Print the dose.
 */
void StepImpl::printDose(PrioritizedMain::Vector *w) {

	int numTarget = mtlb->getNumTargets();
	PrioritizedMain::Vector **D = new PrioritizedMain::Vector *[numTarget];
	//Matrices *D = new Matrices(mtlb->getTTLTargetVoxel(), 1);
	this->getTargetDose(w, D);
	//	int start = 0;
	//	int end = 0;
	int id = 0, i = 0;
	double v;
	printf("Print Dose is called\n");
	int *numvoxel = mtlb->getNumVoxel();
	for (i = 0; i < mtlb->getNumTargets(); i++) {
		//printf("Dose of Target:%d\n",id);
		//D->printNonZeroEntry();
		id = mtlb->getTargetsId(i);
		//end = start + numvoxel[id];
		v = D[i]->getMin();
		printf("Minimum Dose:%lf\n", v);
		v = D[i]->getMax();
		printf("Maximum Dose:%lf\n", v);
		v = D[i]->average();
		printf("Average Dose:%lf\n", v);
		//start = end;
	}
	//	for (int i = 0; i < numTarget; i++) {
	//		delete D[i];
	//	}
	//	delete[] D;

}

double StepImpl::getMinDose(int tar) {

	if (minDose[tar] >= 0)
		return minDose[tar];
	PrioritizedMain::Vector **Dose;
	int numTarget = mtlb->getNumTargets();

	Dose = new PrioritizedMain::Vector *[numTarget];//Matrices(mtlb->getTTLTargetVoxel(), 1);
	getTargetDose(w1, Dose);

	for (int i = 0; i < numTarget; i++) {

		//printf("minDoses\n");
		minDose[i] = Dose[i]->getMin();
		//printf("%lf\n", minDose[i]);

	}
	for (int i = 0; i < numTarget; i++) {
		delete Dose[i];
	}
	delete[] Dose;
	//delete w1;

	return minDose[tar];
}
double StepImpl::getMaxDose(int tar) {
	if (maxDose[tar] >= 0)
		return maxDose[tar];
	PrioritizedMain::Vector **Dose;
	int numTarget = mtlb->getNumTargets();
	Dose = new PrioritizedMain::Vector *[numTarget];//Matrices(mtlb->getTTLTargetVoxel(), 1);
	getTargetDose(w1, Dose);
	for (int i = 0; i < numTarget; i++) {
		//printf("maxDoses\n");
		maxDose[i] = Dose[i]->getMax();
		//printf("%lf\n", maxDose[i]);
	}
	for (int i = 0; i < numTarget; i++) {
		delete Dose[i];
	}
	delete[] Dose;
	return maxDose[tar];
}

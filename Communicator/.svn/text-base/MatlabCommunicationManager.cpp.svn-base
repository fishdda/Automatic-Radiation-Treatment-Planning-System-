/*
 * MatlabCommunicationManager.cpp
 * This class is used for reading Influence Matrices and other configuration data from Matlab planC structure.
 * We used Matlab engine to communicate with the matlab. This class requires that matlab is installed in the system.
 *
 *  Created on: Nov 7, 2011
 *      Author: Paras Babu Tiwari
 */

#define HASHMAPSIZE 4
#include "MatlabCommunicationManager.h"

MatlabCommunicationManager::MatlabCommunicationManager(int caseNum) {
	// TODO Auto-generated constructor stub
	ep = NULL;
	targets = NULL;
	stepIOar = NULL;
	stepIIOar = NULL;
	bndryVoxels = NULL;
	mxStepIIIOar = NULL;
	stepIIorganId = NULL;
	totalTargetRow = 0;
	alpha = NULL;
	//prescDose = 70;
	maxPrescDose = 0;
	ubWeight = Constant::MAXBMLTWT;
	mxHessRow = mxHessCol = mxHessVal = mxRowSize = NULL;
	subMatrix=NULL;
	init(caseNum);

}
//double ** MatlabCommunicationManager::createDoubleMatrix(unsigned int rows,
//		unsigned int cols) {
//	double ** matrix;
//	matrix = (double **) calloc(rows, sizeof(double *));
//	for (unsigned int i = 0; i < rows; i++)
//		matrix[i] = (double *) calloc(cols, sizeof(double));
//	return matrix;
//}
MatlabCommunicationManager::~MatlabCommunicationManager() {
	// TODO Auto-generated destructor stub

	delete numVoxels;
	//printf("Deleted num Voxels Row\n");
//	for (int i = 0; i < maxOrganId; i++) {
//		if (voxels[i] != NULL)
//			delete voxels[i];
//	}

//printf("deleted voxles' content\n");
	delete voxels;
	//printf("deleted voxles\n");
//	for (int i = 0; i < numStepIIOar; i++) {
//		if (alpha[i] != NULL)
//			delete Ä[i];
//	}
	delete alpha;
	//delete targetsId;
	//printf("deleted targets id\n");
	//delete influence;Ä
	//printf("deleted influence matrix\n");
	//delete stepIorganId;
	//printf("deleted stepIOrganId \n");
	if (stepIIorganId == NULL) {
		delete stepIIorganId;
		stepIIorganId = NULL;
	}
	if (alpha == NULL) {
		delete alpha;
		alpha = NULL;
	}

}
/*
 * Get the dose associate with the ids.
 * @param:
 * 	D: Output Parameter, holds the dose value.
 * 	x: Beamlets.
 * 	StrId: Id of the organs.
 * 	numOrgans: Number of organs.
 * 	@return: The dose matrix in the D.
 */
//void MatlabCommunicationManager::getDose(Matrices *D, const double *x,
//		int *StrId, int numOrgans) {
void MatlabCommunicationManager::getDose(PrioritizedMain::Vector ** D,
		const double *x, int *StrId, int numOrgans) {

	assert(D != NULL);
	int offset = 0;
	PrioritizedMain::Vector *w = new PrioritizedMain::Vector(x,
			influence->getCol());
	getDose(D, w, StrId, numOrgans);
	delete w;
}
/*
 * Get the dose associate with the ids.
 * @param:
 * 	D: Output Parameter, holds the dose value.
 * 	w: Beamlets.
 * 	StrId: Id of the organs.
 * 	numOrgans: Number of organs.
 * 	@return: The dose matrix in the D.
 */

 void MatlabCommunicationManager::getDose(PrioritizedMain::Vector **D,
		PrioritizedMain::Vector *w, int *StrId, int numOrgans) {

	assert(D != NULL);
	int offset = 0;
	
	//printf("numOrgans:%d\n",numOrgans);
	for (int i = 0; i < numOrgans; i++) {

		int id = StrId[i];
		//printf("id:%d\n",id);
		if(subMatrix==NULL)
		{
			printf("Inside submatrix\n");
			subMatrix = new Matrices*[numOrgans];
			subMatrix[i] = influence->getSubMatrix(id, voxels[id], numVoxels[id]);
		//Matrices *result = new Matrices(subMatrix->getRow(), w->getCol());
		
		}
		D[i] = subMatrix[i]->matrixMultiply(w);
		
		//delete result;

	}
}

//void MatlabCommunicationManager::getDose(PrioritizedMain::Vector **D,
//		PrioritizedMain::Vector *w, int *StrId, int numOrgans) {
//
//	assert(D != NULL);
//	int offset = 0;
//	Matrices *subMatrix;
//	//printf("numOrgans:%d\n",numOrgans);
//	for (int i = 0; i < numOrgans; i++) {
//
//		int id = StrId[i];
//		//printf("id:%d\n",id);
//		subMatrix = influence->getSubMatrix(id, voxels[id], numVoxels[id]);
//		//Matrices *result = new Matrices(subMatrix->getRow(), w->getCol());
//		D[i] = subMatrix->matrixMultiply(w);
//		delete subMatrix;
//		//delete result;
//
//	}
//}
void MatlabCommunicationManager::init(int caseNum) {
	int ret;
	if (!(ep = engOpen("\0"))) {
		fprintf(stderr, "\nCan't start MATLAB engine\n");

	} else {

		putCaseId(caseNum);
		ret = engEvalString(ep, "getInfluenceMatrix");

	}
	caseNo = caseNum;
	//printf("Calling getConfigParam\n");
	getConfigParam();
	//printf("Calling makespace\n");
	makeSpace();
	//printf("Calling getTargetVoxel\n");
	getTargetVoxel();
	//printf("Calling getInfluenceMatrix\n");

	getInfluenceMatrix();
	//printf("Calling getStepIOarVoxels\n");
	getStepIOarVoxels();
	//printf("Calling computerNonZeroVoxelofTargets\n");
	computerNonZeroVoxelofTargets();
	//printf("Calling getHessianMatrix\n");
	getHessianMatrix();
	//printf("Calling getStepIIVoxels\n");
	getStepIIVoxels();
	hessians = (SparseMatrix **) calloc(numTargets, sizeof(SparseMatrix *));
	getHessians();
	getStepIIIInfo();

	//freeDblArray(bndryVoxels);
	//freeDblArray(mxHessCol);
	//freeDblArray(mxHessRow);
	freeDblArray(mxHessVal);
	engClose(ep);
	//analyzeAccess();
//	printf("# of targets:%d\n", numTargets);
//	printf("# of stepI organ:%d\n", numStepIOar);
//	printf("# step II organs:%d\n", numStepIIOar);
//	printf("Target id\n");
//	for (int i = 0; i < numTargets; i++)
//		printf("%lf  ", targetsId[i]);
//
//	printf("\nStepI organ Id\n");
//	for (int i = 0; i < numStepIOar; i++)
//		printf("%lf ", stepIorganId[i]);
//	printf("\nStepII organ id\n");
//	for (int i = 0; i < numStepIIOar; i++)
//		printf("%d ", getStepIId(i));
//	printf("\nAlpha\n");
//	for (int i = 0; i < numStepIIOar; i++)
//		printf("%lf ", alpha[i]);
//	printf("# of voxels\n");
//	for (int i = 0; i < maxOrganId; i++)
//		printf("%d  ", numVoxels[i]);
}
void MatlabCommunicationManager::breakInfMatrix()
{
	
	targetInf = new Matrices*[numTargets];
	for (int i = 0; i < numTargets; i++) {
		int id = targetsId[i];
		targetInf[i] = influence->getSubMatrix(id,voxels[id],
				numVoxels[id]);
	}
	stepIinf = new Matrices*[numStepIOar];
	for (int i = 0; i < numStepIOar; i++) {
		int id = stepIorganId[i];
		stepIinf[i] = influence->getSubMatrix(id, voxels[id],
				numVoxels[id]);
	}
}
void MatlabCommunicationManager::analyzeAccess() {

	time_t start_t, end_t;
	start_t = time(NULL);
	double d;
	for (long unsigned int i = 0; i < influence->getRow(); i++) {
		for (long unsigned int j = 0; j < influence->getCol(); j++) {
			d = influence->getAt(i, j);
		}
	}
	end_t = time(NULL);
	double diff = difftime(end_t, start_t);
	printf("Matlab Row access elapsed time:%lf seconds\n", diff);
	start_t = time(NULL);
	for (long unsigned int j = 0; j < influence->getCol(); j++) {
		for (long unsigned int i = 0; i < influence->getRow(); i++) {

			d = influence->getAt(i, j);
		}
	}
	end_t = time(NULL);
	diff = difftime(end_t, start_t);
	printf("Matlab column access elapsed time:%lf seconds\n", diff);
}
void MatlabCommunicationManager::getStepIIVoxels() {
	//printf("Inside stepII voxel\n");
	if (stepIIOar == NULL) {
		printf("Getting stepII voxel\n");
		stepIIOar = engGetVariable(ep, "step2OARs");
		if (stepIIOar == NULL)
			printf("Sorry You didn't create step2OARs.\n");
	}
	double *organId;
	//printf("Got stepII voxel\n");
	getOarVoxels(stepIIOar, &organId, numStepIIOar, ttlStepIIOarRow);
	//printf("Number of StepII Oar:%d\n", numStepIIOar);
	//printf("Total # of stepII organ's voxels:%d\n", ttlStepIIOarRow);
	if (stepIIorganId == NULL)
		stepIIorganId = new int[numStepIIOar];
	for (int i = 0; i < numStepIIOar; i++) {
		stepIIorganId[i] = (int) organId[i];
		//printf("stepIIorganId[%d]:%d\n", i, stepIIorganId[i]);
	}
	//delete[] organId;
	if ((mxAlpha = engGetVariable(ep, "alpha")) == NULL) {
		printf("Oops! You didn't create a variable step2XVals.\n\n");
	}
	//alpha = getDoubleArray(mxAlpha);
	if (alpha == NULL) {

		alpha = new double*[numStepIIOar];
	}
	const mxArray *cell_element_ptr;
	for (int i = 0; i < numStepIIOar; i++) {
		cell_element_ptr = mxGetCell(mxAlpha, i);
		alpha[i] = getDoubleArray(cell_element_ptr);
	}

	for (int i = 0; i < numStepIIOar; i++) {

		printf("alpha:%lf\n", alpha[i][0]);
	}
	getNzStepIIVoxl();
}
/*
 * Get the stepIII organs id and the voxels corresponding to these organs.
 * @Author: Paras Babu Tiwari
 */
void MatlabCommunicationManager::getStepIIIInfo() {
	if (mxStepIIIOar == NULL) {
		mxStepIIIOar = engGetVariable(ep, "step3OARs");
		if (mxStepIIIOar == NULL)
			printf("Sorry You didn't create step3OARs.\n");
	}
	double *organId;

	getOarVoxels(mxStepIIIOar, &organId, numStepIIIOar, ttlStepIIIOarRow);
	//printf("Number of StepIII Oar:%d\n", numStepIIIOar);
	//printf("Total # of stepIII organ's voxels:%d\n", ttlStepIIIOarRow);
	stepIIIorganId = new int[numStepIIIOar];
	for (int i = 0; i < numStepIIIOar; i++) {
		stepIIIorganId[i] = (int) organId[i];
		//printf("stepIIIorganId[%d]:%d\n", i, stepIIIorganId[i]);
	}
	//delete[] organId;

	//getNzStepIIVoxl();
}
/*
 * Get the # of non-zero voxel of the stepII organs. Modifies class level variable.
 * @Input:None
 * @return: None
 * @Author: Paras Babu Tiwari
 */
void MatlabCommunicationManager::getNzStepIIVoxl() {

	const mxArray* nonzero;
	nonzero = engGetVariable(ep, "numStepII");
	if (nonzero != NULL) {
		double *ptr = (double *) mxGetPr(nonzero);
		nonzeroStepIIVoxels = (int) *ptr;
		//printf("nonzerostepIVx:%d", nonzeroStepIVx);
		//delete ptr;
	}
}
void MatlabCommunicationManager::getConfigParam() {

	mxArray *morganid;
	if ((morganid = engGetVariable(ep, "morganId")) == NULL) {
		//printf("morganid:%x\n", morganid);
		printf("Oops! You didn't create a variable morganid.\n\n");
	}
	maxOrganId = (int) getDoubleVal(morganid);

	maxOrganId = maxOrganId + 1;
	freeDblArray(morganid);
	printf("Max Organ Id:%d\n", maxOrganId);
}
void MatlabCommunicationManager::makeSpace() {

	numVoxels = new int[maxOrganId];

	voxels = new double*[maxOrganId];
	for (int i = 0; i < maxOrganId; i++)
		voxels[i] = NULL;
}
void MatlabCommunicationManager::computerNonZeroVoxelofTargets() {

	double *dblpr;
	int value = 0;
	mxArray *nonZeroVoxel = NULL;
	if (nonZeroVoxel == NULL) {
		nonZeroVoxel = engGetVariable(ep, "num");

	}
	//printf("Nonzero voxel:%x",nonZeroVoxel);
	if (nonZeroVoxel != NULL)
		dblpr = (double *) mxGetPr(nonZeroVoxel);

	if (dblpr != NULL)
		value = (int) (*dblpr);
	nonzeroTargetVx = value;
	freeDblArray(nonZeroVoxel);
}

long unsigned int MatlabCommunicationManager::getRowInfMat() {
	return influence->getRow();
}
int MatlabCommunicationManager::getColInfMat() {

	return influence->getCol();
}
double *MatlabCommunicationManager::getPrescribedDose() {

	return prescDose;

}
double MatlabCommunicationManager::getMaxPrescDose() {
	double max = prescDose[0];
	for (int i = 0; i < numTargets; i++) {
		if (max < prescDose[i])
			max = prescDose[i];
	}
	return max;
}
double MatlabCommunicationManager::getMaxDose() {
	if (maxPrescDose == 0)
		maxPrescDose = getMaxPrescDose();
	//double max = 1.15 * maxPrescDose;
	double max = 16;
	//printf("max dose:%lf\n",max);
	return max;
}
double MatlabCommunicationManager::getUBWeight() {

	return ubWeight;
}

double MatlabCommunicationManager::getDoubleVal(const mxArray *ptr) {
	double *dblptr = NULL;
	double val;
	if (ptr != NULL)
		dblptr = (double *) mxGetPr(ptr);
	if (dblptr != NULL)
		val = (double) (*dblptr);
	//delete dblptr;
	return val;

}
double *MatlabCommunicationManager::getDoubleArray(const mxArray *ptr) {
	double *dblptr = NULL;
	if (ptr != NULL)
		dblptr = (double *) mxGetPr(ptr);
	return dblptr;
}
void MatlabCommunicationManager::freeDblArray(mxArray *ptr) {
	mxDestroyArray(ptr);
}
double ** MatlabCommunicationManager::getDoubleArrays(const mxArray *ptr) {
	double **dblptr = NULL;
	if (ptr != NULL)
		dblptr = (double **) mxGetPr(ptr);
	return dblptr;
}
int MatlabCommunicationManager::getNumTargets() {
	int total_num_of_elements = 0;
	if (targets == NULL) {
		if (ep != NULL) {
			targets = engGetVariable(ep, "targets");
		}
	}
	if (targets != NULL) {
		total_num_of_elements = mxGetNumberOfElements(targets);
	}

	return total_num_of_elements;
}

/*
 * This function is used to get the hessian matrix of each targets
 */
void MatlabCommunicationManager::getHessians() {
	mxArray *rowCell, *colCell, *valCell;
	double *rowSize;
	double *row, *col, *val;
	if ((mxRowSize = engGetVariable(ep, "rowSize")) == NULL) {
		printf("Oops! You didn't create a variable rowsize.\n\n");

	}
	rowSize = getDoubleArray(mxRowSize);
	if (mxHessRow == NULL) {
		mxHessRow = engGetVariable(ep, "hessRow");
		if (mxHessRow == NULL)
			printf("Can't read hessRow\n");
	}
	if (mxHessCol == NULL) {
		mxHessCol = engGetVariable(ep, "hessCol");
		if (mxHessCol == NULL)
			printf("Can't read hessCol\n");
	}
	if (mxHessVal == NULL) {
		mxHessVal = engGetVariable(ep, "hessVal");
		if (mxHessVal == NULL)
			printf("Can't read hessVal\n");
	}

	for (int i = 0; i < numTargets; i++) {
		printf("row size:%lf\n", rowSize[i]);
	}
	int n;
	for (int i = 0; i < numTargets; i++) {
		rowCell = mxGetCell(mxHessRow, i);

		row = getDoubleArray(rowCell);
		n = mxGetNumberOfElements(rowCell);
		//printf("# of row:%d\n", n);
//		for (int j = 0; j < 10; j++)
//			printf("row[%d]:%lf\n", j, row[j]);

		colCell = mxGetCell(mxHessCol, i);
		n = mxGetNumberOfElements(colCell);
		printf("# of col:%d\n", n);
		col = getDoubleArray(colCell);
//		for (int j = 0; j < 10; j++)
//			printf("col[%d]:%lf\n", j, col[j]);

		valCell = mxGetCell(mxHessVal, i);
		val = getDoubleArray(valCell);
		n = mxGetNumberOfElements(colCell);
		//printf("# of val:%d\n", n);

		hessians[i] = new SparseMatrix(row, col, val, rowSize[i],
				this->getColInfMat(), this->getColInfMat());
		//
		//printf("hessians:%x\n",hessians[i]);
	}
	//freeDblArray(mxHessVal);
}
/**
 * Get the voxels of Organs' id given in oar.
 * Input Param:
 * 	@oar: The organ id in mx format. This should be non-null.
 * Output Param:
 * 	@oarId: Convert the mx format oar to double array.
 * 	@numOar: Number of organs in this step;
 * 	@ttlRow: Total number of row in this step.
 */
void MatlabCommunicationManager::getOarVoxels(mxArray *oar, double **oarId,
		int& numOar, int& ttlRow) {
	int i = 0;

	const mxArray *cell_element_ptr;
	const mxArray* nonzero;
	double *pr, *inpr;
	int index;
	if (bndryVoxels == NULL) {
		bndryVoxels = engGetVariable(ep, "allVoxelC");
	}
	if (oar == NULL) {
		printf("Missing organs id.\n");
		return;
	}
	//printf("oarid Address:%x\n", *oarId);
	*oarId = getDoubleArray(oar);
	//printf("oarid Address:%x\n", *oarId);
	//printf("stepIIorganId Address:%x\n", stepIIorganId);
	numOar = (int) mxGetNumberOfElements(oar);
	inpr = (double *) mxGetPr(oar);
	double *orig = inpr;
//	printf("Number of Step 2 oar:%d\n", numOar);
	for (i = 1; i <= numOar; i++) {
		printf("id:%lf\n", (*oarId)[i - 1]);
		index = (int) (*inpr);
		index = index - 1;
		cell_element_ptr = mxGetCell(bndryVoxels, index);
		voxels[index + 1] = getDoubleArray(cell_element_ptr);
		//printf("index:%d\n",index);
		numVoxels[index + 1] = (int) mxGetNumberOfElements(cell_element_ptr);
		printf("Voxels of organ:%d:%d\n", index + 1, numVoxels[index + 1]);
		//double *vx = getDoubleArray(cell_element_ptr);
		//voxels[index + 1] = new double[numVoxels[index+1]];
		//printf("vox Address:%x\n",voxels[index+1]);
		//double *r = voxels[index+1];
		//double *e = vx+ numVoxels[index+1];
		//std::copy(vx,e,r);

		ttlRow += (int) mxGetNumberOfElements(cell_element_ptr);
		inpr++;
		//totVoxels += mxGetNumberOfElement(mxGetData(bndryVoxels,index));
	}
//	printf("numStepIIOar:%d\n",numStepIIOar);
//	for (int i = 0; i < numStepIIOar; i++) {
//		printf("stepIIorganId[%d]:%lf", i, stepIIorganId[i]);
//	}
//	printf("Total Row:%d\n", ttlRow);
	inpr = orig;

//	nonzero = engGetVariable(ep, "numStepI");
//	if (nonzero != NULL) {
//		double *ptr = (double *) mxGetPr(nonzero);
//		nonzeroStepIVx = (int) *ptr;
//		printf("nonzerostepIVx:%d", nonzeroStepIVx);
//		//delete ptr;
//	}
//	//delete inpr;
}

void MatlabCommunicationManager::getStepIOarVoxels() {
	int i = 0;
	numStepIOar = 0;
	ttlStepIOarVoxel = 0;
	const mxArray *cell_element_ptr;
	mxArray* nonzero;
	double *pr, *inpr;
	int index;
	if (bndryVoxels == NULL) {
		bndryVoxels = engGetVariable(ep, "allVoxelC");
	}
	if (stepIOar == NULL) {
		stepIOar = engGetVariable(ep, "step1OARs");
	}
	//stepIorganId = getDoubleArray(stepIOar);
	double *id = getDoubleArray(stepIOar);

	numStepIOar = (int) mxGetNumberOfElements(stepIOar);
	inpr = (double *) mxGetPr(stepIOar);
	double *orig = inpr;
	stepIorganId = new int[numStepIOar];
	for (i = 1; i <= numStepIOar; i++) {

		index = (int) (*inpr);
		index = index - 1;
		cell_element_ptr = mxGetCell(bndryVoxels, index);
		voxels[index + 1] = getDoubleArray(cell_element_ptr);

		//printf("index:%d\n",index);
		numVoxels[index + 1] = (int) mxGetNumberOfElements(cell_element_ptr);
		printf("Voxels of Step I organ:%d:%d\n", index + 1,
				numVoxels[index + 1]);

		ttlStepIOarVoxel += (int) mxGetNumberOfElements(cell_element_ptr);
		inpr++;
		stepIorganId[i - 1] = id[i - 1];
		//totVoxels += mxGetNumberOfElement(mxGetData(bndryVoxels,index));
	}
	inpr = orig;
	nonzero = engGetVariable(ep, "numStepI");
	if (nonzero != NULL) {
		double *ptr = (double *) mxGetPr(nonzero);
		nonzeroStepIVox = (int) *ptr;
		printf("nonzerostepIVx:%d", nonzeroStepIVox);
		//delete ptr;
	}
	//freeDblArray(cell_element_ptr);
	freeDblArray(nonzero);
//	printf("nonzerostepIVx:%d", nonzeroStepIVox);
	//delete inpr;
}

void MatlabCommunicationManager::getTargetVoxel() {
	const mxArray *cell_element_ptr;
	int index;

	if (targets == NULL) {
		if (ep != NULL) {
			targets = engGetVariable(ep, "targets");
		}
	}
	targetsId = getDoubleArray(targets);
	numTargets = (int) mxGetNumberOfElements(targets);
	if (bndryVoxels == NULL) {
		bndryVoxels = engGetVariable(ep, "allVoxelC");
	}
	//printf("tid:%lf\n",targetsId[0]);
	double *tId = targetsId;
	//printf("Number of Targets:%d\n", numTargets);
	for (int i = 0; i < numTargets; i++) {
		//printf("Inside the loop\n");
//		printf("tId:%x\n",tId);
		index = (int) (*targetsId);
		index = index - 1;
		//printf("Index:%d\n", index);
		cell_element_ptr = mxGetCell(bndryVoxels, index);
		int n = (int) mxGetNumberOfElements(cell_element_ptr);
		//printf("index:%d,n:%d\n", index + 1, n);
		numVoxels[index + 1] = n;
		totalTargetRow += n;
		voxels[index + 1] = getDoubleArray(cell_element_ptr);

		printf("Voxels of target:%d:%d\n", index + 1, n);

		//printf("voxels:%x\n", voxels[index + 1]);
		double *row = voxels[index + 1];
		//printf("row:%x\n", row);
		//printf("\n");
		targetsId++;

	}

	targetsId = tId;
	getPrescDose();
	//printf("tid:%lf\n",targetsId[0]);

}
void MatlabCommunicationManager::getPrescDose() {
	const mxArray *p;
	p = engGetVariable(ep, "Dpre");
	if (p == NULL) {
		printf("Oops! You didn't create a variable dosePrescribedTarget.\n\n");
	}
	prescDose = getDoubleArray(p);

}
void MatlabCommunicationManager::analyze(const mxArray *array_ptr) {
	mxClassID category;

	category = mxGetClassID(array_ptr);

	switch (category) {
	case mxLOGICAL_CLASS:
		printf("Logical\n");
		break;
	case mxCHAR_CLASS:
		printf("Char\n");
		break;
	case mxSTRUCT_CLASS:
		printf("Struct\n");
		break;
	case mxCELL_CLASS:
		printf("Cell\n");
		break;
	case mxUNKNOWN_CLASS:
		printf("Unknown Class");
		break;
	default:
		analyze_full(array_ptr);
		break;

	}

}

void MatlabCommunicationManager::analyze_full(
		const mxArray *numeric_array_ptr) {
	mxClassID category;

	category = mxGetClassID(numeric_array_ptr);
	switch (category) {
	case mxINT8_CLASS:
		printf("Int8\n");
		break;
	case mxUINT8_CLASS:
		printf("Uint8\n");
		break;
	case mxINT16_CLASS:
		printf("int16\n");
		break;
	case mxUINT16_CLASS:
		printf("uint16\n");
		break;
	case mxINT32_CLASS:
		printf("int32\n");
		break;
	case mxUINT32_CLASS:
		printf("uint32\n");
		break;
	case mxINT64_CLASS:
		printf("int64\n");
		break;
	case mxUINT64_CLASS:
		printf("uint64\n");
		break;
	case mxSINGLE_CLASS:
		printf("single\n");
		break;
	case mxDOUBLE_CLASS:
		printf("Double\n");
		break;
	default:
		break;
	}
}
void MatlabCommunicationManager::getInfluenceMatrix() {

	mxArray *r, *c, *v, *nrow;
	//printf("1\n");
	//

	printf("ep:%x\n", ep);
	if ((r = engGetVariable(ep, "r")) == NULL) {
		printf("r:%x\n", r);
		printf("Oops! You didn't create a variable r.\n\n");

	}

	if ((c = engGetVariable(ep, "c")) == NULL) {
		printf("Oops! You didn't create a variable c.\n\n");
	}

	if ((v = engGetVariable(ep, "val")) == NULL) {
		printf("Oops! You didn't create a variable val.\n\n");
	}

	if ((nrow = engGetVariable(ep, "numRow")) == NULL) {
		printf("Oops! You didn't create a variable numRow.\n\n");
	}
	//analyze(r);
	//printf("2\n");
	double *nr = getDoubleArray(r);
	//printf("1.1\n");
	double *nc = getDoubleArray(c);
	//printf("1.2\n");
	double *val = getDoubleArray(v);
	//printf("1.3\n");
	long unsigned int numRow = (long unsigned int) getDoubleVal(nrow);
	//printf("Number of row:%lf\n",numRow);
	//printf("Number of rows and col:%lf %lf\n",nr,nc);
//	printf("printing");
//	printf("%lf\n",nr[0]);
	const mxArray *rf, *cf;
	if ((rf = engGetVariable(ep, "row")) == NULL) {
		printf("Oops! You didn't create a variable row.\n\n");

	}
	if ((cf = engGetVariable(ep, "col")) == NULL) {
		printf("Oops! You didn't create a variable col.\n\n");

	}

	long unsigned int ncf, nrf;

	ncf = (long unsigned int) getDoubleVal(cf);
	nrf = (long unsigned int) getDoubleVal(rf);
	printf("Number of rows and col:%lu %lu\n", nrf, ncf);
	//int n = getNumStepIOar() + getNumTargets();
	//Because the Hessian has the term for target also.

	influence = new SparseMatrix(nr, nc, val, numRow, nrf, ncf); //Need to fix the constant value
	//mxArray *nvar;
//	if ((nvar = engGetVariable(ep, "numVar")) == NULL) {
//		printf("Oops! You didn't create a variable c.\n\n");
//	}
	//printf("Got non zero beamlete\n");
	//nzbmltfortar = (int) getDoubleVal(nvar);
	//printf("Non zero beamlet:%d\n",nzbmltfortar);
	freeDblArray(v);
	//freeDblArray(r);
	//freeDblArray(c);

}
/*
 * Get the voxels corresponding to the organ ids provided in the oarId.
 * @param:
 * 	oarId: The organ  id.
 * 	numOrgan: The number of elements in the organ.
 * 	@return the
 */
void MatlabCommunicationManager::getVoxelofOar(int *oarId, int numOrgan,
		double *vxls) {
	int k = 0;
	for (int i = 0; i < numOrgan; i++) {
		for (int j = 0; j < numVoxels[oarId[i]]; j++) {
			vxls[k] = voxels[oarId[i]][j];
			k++;
		}
	}
}
/*
 * Get the Hessian Matrices. The matrix indicies are stored in the matlab format. It's caller responsibility to convert in C format.
 */

void MatlabCommunicationManager::getHessianMatrix() {

	mxArray *r, *c, *v, *nrow;
	//printf("1\n");

	if ((r = engGetVariable(ep, "lr")) == NULL) {
		printf("Oops! You didn't create a variable lr.\n\n");

	}

	if ((c = engGetVariable(ep, "lc")) == NULL) {
		printf("Oops! You didn't create a variable lc.\n\n");
	}

	if ((v = engGetVariable(ep, "lv")) == NULL) {
		printf("Oops! You didn't create a variable lv.\n\n");
	}

	//analyze(r);
	//printf("2\n");
	double *nr = getDoubleArray(r);
	//printf("1.1\n");
	double *nc = getDoubleArray(c);
	//printf("1.2\n");
	double *val = getDoubleArray(v);
	//printf("1.3\n");

	//printf("Number of row:%lf\n",numRow);
	//printf("Number of rows and col:%lf %lf\n",nr,nc);
//	printf("printing");
//	printf("%lf\n",nr[0]);
	const mxArray *rf, *cf;
	if ((nrow = engGetVariable(ep, "nlr")) == NULL) {
		printf("Oops! You didn't create a variable r.\n\n");

	} else {
		printf("Got the variable nlr\n");
	}
	long unsigned int numRow = (long unsigned int) getDoubleVal(nrow);
	printf("Number of rows:%lu\n", numRow);
	int n = this->getColInfMat() + getNumTargets();

	hessian = new SparseMatrix(nr, nc, val, numRow, n, n);
	freeDblArray(v);
	//freeDblArray(r);
	//freeDblArray(c);

}
void MatlabCommunicationManager::putCaseId(int caseNum) {
	const mxArray *pm;
	char *name = "patID";
	pm = mxCreateDoubleScalar(caseNum);
	engPutVariable(ep, name, pm);
}
double *MatlabCommunicationManager::getMetrics(PrioritizedMain::Vector *w1,
		PrioritizedMain::Vector *w2, PrioritizedMain::Vector *w3,
		PrioritizedMain::Vector *w4) {
	const mxArray *pm1, *mxD95, *pm2, *pm3, *pm4;
	double *data1, *data2, *data3, *data4;
	const char *name1 = "wVI";
	const char *name2 = "wVII";
	const char *name3 = "wVIII";
	const char *name4 = "wVIV";
	pm1 = mxCreateDoubleMatrix(w1->size(), 1, mxREAL);
	pm2 = mxCreateDoubleMatrix(w2->size(), 1, mxREAL);
	pm3 = mxCreateDoubleMatrix(w3->size(), 1, mxREAL);
	pm4 = mxCreateDoubleMatrix(w4->size(), 1, mxREAL);
	if (pm1 == NULL) {
		printf("Can't create an array\n");
		exit(0);
	}
	data1 = mxGetPr(pm1);
	data2 = mxGetPr(pm2);
	data3 = mxGetPr(pm3);
	data4 = mxGetPr(pm4);
	getData(data1, w1);
	getData(data2, w2);
	getData(data3, w3);
	getData(data4, w4);

	printf("Got data\n");
	if (!(ep = engOpen("\0"))) {
		fprintf(stderr, "\nCan't start MATLAB engine\n");

	} else {

		putCaseId(caseNo);

		engPutVariable(ep, name1, pm1);
		engPutVariable(ep, name2, pm2);
		engPutVariable(ep, name3, pm3);
		engPutVariable(ep, name4, pm4);
		engEvalString(ep, "calcD95");
	}
	if ((mxD95 = engGetVariable(ep, "metrics")) == NULL) {
		printf("Oops! You didn't create a variable metrics.\n\n");

	}
	engClose(ep);
	double *metrics = getDoubleArray(mxD95);
	//printf("d95[0]:%lf\n", d95[0]);
	return metrics;
}
void MatlabCommunicationManager::getData(double *d,
		PrioritizedMain::Vector *w) {
	for (int i = 0; i < w->size(); i++) {
		d[i] = w->getAt(i);
	}
}

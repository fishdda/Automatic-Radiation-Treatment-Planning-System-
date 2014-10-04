/*
 * oneVoxel.cpp
 *
 *  Created on: Dec 28, 2011
 *      Author: user
 */

#include "oneVoxel.h"

oneVoxel::oneVoxel() {
	// TODO Auto-generated constructor stub
	printf("Super executed\n");

	voxels = (double **) calloc(13, sizeof(double *));

	voxels[2] = (double *) calloc(1, sizeof(double));
	voxels[10] = (double *) calloc(1, sizeof(double));
	voxels[11] = (double *) calloc(1, sizeof(double));
	voxels[12] = (double *) calloc(1, sizeof(double));
	stepIorganId = new int[1];
	targetsId = new int[3];
	stepIorganId[0] = 1;
	targetsId[0] = 10;
	targetsId[1] = 11;
	targetsId[2] = 12;

	numVoxelsRow = new int[13];
	numStepIOar = 1;
	numStepIIIOar = 1;
	stepIIorganId = new int[1];
	stepIIorganId[0] = 2;
	stepIIIorganId = new int[1];
	stepIIIorganId[0] = 3;
	//stepIIIorganId[0] = 4;
	init();
	prescDose = new double[numTargets];
	for(int i=0;i<numTargets;i++)
		prescDose[i]= 70;

}

oneVoxel::~oneVoxel() {
	// TODO Auto-generated destructor stub
}
/*
 * Get the dose associate with the ids.
 * @param:
 * 	D: Output Parameter, holds the dose value. Should be not null.
 * 	x: Beamlets.
 * 	StrId: Id of the organs.
 * 	numOrgans: Number of organs.
 * 	@return: The dose matrix in the D.
 */
void oneVoxel::getDose(PrioritizedMain::Vector ** D, const double *x, int *StrId, int numOrgans) {

	assert(D != NULL);
	PrioritizedMain::Vector *w = new PrioritizedMain::Vector(x, influence->getCol());
	getDose(D,w,StrId,numOrgans);
	delete w;
//	for (int i = 0; i < numOrgans; i++) {
//
//		int id = StrId[i];
//		subMatrix = influence->getSubMatrix(id, voxels[id], numVoxelsRow[id]);
//		//printf("SubMatrix\n");
//		//subMatrix->printNonZeroEntry();
//		Matrices *result = new Matrices(subMatrix->getRow(), w->getCol());
//		result->matrixMultiply(subMatrix, w);
//		D->copy(result, offset);
//		offset += getNumVoxels(id);
//		delete result;
//
//	}
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
void oneVoxel::getDose(PrioritizedMain::Vector **D, PrioritizedMain::Vector *w, int *StrId, int numOrgans) {

	assert(D != NULL);
	int offset = 0;
	Matrices *subMatrix;

	for (int i = 0; i < numOrgans; i++) {

		int id = StrId[i];
		subMatrix = influence->getSubMatrix(id, voxels[id], getNumVoxels(id));
		//printf("submatrix of %d\n",id);
		//subMatrix->printNonZeroEntry();
		//Matrices *result = new Matrices(subMatrix->getRow(), w->getCol());
		D[i] = subMatrix->matrixMultiply(w);
		//delete result;

	}
}
void oneVoxel::init() {

	getTargetVoxel();
	getInfluenceMatrix();
	getVoxelStepI();
	ttlTargetRow = 0;
	for (int i = 0; i < numTargets; i++) {
		ttlTargetRow += numVoxelsRow[targetsId[i]];
	}
	ttlStepIOarRow = 0;
	for (int i = 0; i < numStepIOar; i++) {
		ttlStepIOarRow += numVoxelsRow[stepIorganId[i]];
	}
	getStepIIVoxel();
	printf("before stepiii col:%d\n",influence->getCol());
	getStepIIIinfo();
	printf("after stepiii col:%d\n",influence->getCol());
	//influence->setNumFullRow(500);
	//influence->setNumFullCol(500);
	//printf("numCol:%d\n",numCol);
	//printf("numRow:%d\n",numRow);
}
void oneVoxel::getStepIIIinfo()
{
	double *row;
	row = new double[4];
	for(int i=35;i<39;i++)
	{
		row[i-35] = i;
	}
	numVoxelsRow[stepIIIorganId[0]] = 4;
	voxels[stepIIIorganId[0]] = row;

}
void oneVoxel::getStepIIVoxel() {

	double *row, *row1, *row2;
//	FileManager *f = new FileManager("/Users/user/CoinIpopt/Ipopt/examples/Project/Debug/tarVox.txt");
//	row = (double *)malloc(sizeof(double)*5206);
//	f->readFile(row);
//	delete f;
	//row = (double *) malloc(sizeof(double) * 4);
	row = new double[4];
	for (int i = 5; i < 9; i++) {
		row[i - 5] = i;
	}
	numVoxelsRow[stepIIorganId[0]] = 4;
	voxels[stepIIorganId[0]] = row;
	//alpha = new double[1];
	alpha = new double*[1];
	alpha[0] = new double[1];
	alpha[0][0] = 0.5;

}

void oneVoxel::getTargetVoxel() {

	double *row, *row1, *row2;
//	FileManager *f = new FileManager("/Users/user/CoinIpopt/Ipopt/examples/Project/Debug/tarVox.txt");
//	row = (double *)malloc(sizeof(double)*5206);
//	f->readFile(row);
//	delete f;
	row = (double *) malloc(sizeof(double) * 5);
	row1 = (double *) malloc(sizeof(double) * 10);
	row2 = (double *) malloc(sizeof(double) * 15);
	numTargets = 3;

	for (int i = 0; i < 5; i++) {
		row[i] = i;
	}
	numVoxelsRow[10] = 5;
	voxels[targetsId[0]] = row;

	for (int i = 0; i < 10; i++) {
		row1[i] = i + 10;
	}
	numVoxelsRow[11] = 10;
	voxels[targetsId[1]] = row1;

	for (int i = 0; i < 15; i++) {
		row2[i] = i + 20;
	}
	numVoxelsRow[12] = 15;
	voxels[targetsId[2]] = row2;
	hessians = (SparseMatrix **) calloc(numTargets, sizeof(SparseMatrix *));
}
void oneVoxel::getVoxelStepI() {
	double *row;
	row = (double *) malloc(sizeof(double) * 10);
	for (int i = 39; i < 45; i++) {
		row[i - 39] = i;
	}
//	double *row;
	//row = (double *)malloc(sizeof(double)*671);
//	FileManager *f = new FileManager("/Users/user/CoinIpopt/Ipopt/examples/Project/Debug/stepIVox.txt");
//	f->readFile(row);
//	delete f;
	voxels[stepIorganId[0]] = row;
	numVoxelsRow[stepIorganId[0]] = 6;
}
void oneVoxel::getInfluenceMatrix() {

	double *row, *col, *val, *hessVal, *row1, *col1, *val1, *row2, *col2, *val2,
			*row3, *col3, *val3, *hessRow, *hessCol;
	row = new double[12];
	col = new double[12];
	val = new double[12];
	hessVal = new double[6];
	row1 = new double[2];
	col1 = new double[2];
	val1 = new double[2];
	row2 = new double[2];
	col2 = new double[2];
	val2 = new double[2];
	row3 = new double[2];
	col3 = new double[2];
	val3 = new double[2];
	hessRow = new double[6];
	hessCol = new double[6];
	//printf("Assigning\n");
//	FileManager *f;
//	row = (double *)malloc(sizeof(double)*1262837);
//	f = new FileManager("/Users/user/CoinIpopt/Ipopt/examples/Project/Debug/r.txt");
//	f->readFile(row);
//	delete f;
//
//	col = (double *)malloc(sizeof(double)*1262837);
//	f = new FileManager("/Users/user/CoinIpopt/Ipopt/examples/Project/Debug/c.txt");
//	f->readFile(col);
//	delete f;
//
//	val = (double *)malloc(sizeof(double)*1262837);
//	f = new FileManager("/Users/user/CoinIpopt/Ipopt/examples/Project/Debug/val.txt");
//	f->readFile(val);
//	delete f;

	for (int i =0; i < 2; i++) {
		//totalTargetRow("Assigning\n");
		//printf("row:%d col:%d\n",i,i);
		row[i] = i;
		col[i] = i;
		row1[i] = i;
		col1[i] = i;
		hessRow[i] = i;
		hessCol[i] = i;

	}
	for (int i = 5; i < 7; i++) {
		//totalTargetRow("Assigning\n");
		//printf("row:%d col:%d\n",i,i);
		row[i - 3] = i;
		col[i - 3] = i;

	}
	for (int i = 10; i < 12; i++) {
		//totalTargetRow("Assigning\n");
		//printf("row:%d col:%d\n",i,i);
		row[i - 6] = i;
		col[i - 6] = i;
		row2[i - 10] = i;
		col2[i - 10] = i;
		hessRow[i - 8] = i;
		hessCol[i - 8] = i;
	}
	for (int i = 20; i < 22; i++) {
		//totalTargetRow("Assigning\n");
		//printf("row:%d col:%d\n",i,i);
		row[i - 14] = i;
		col[i - 14] = i;
		row3[i - 20] = i;
		col3[i - 20] = i;
		hessRow[i - 16] = i;
		hessCol[i - 16] = i;

	}
	for (int i = 35; i < 37; i++) {
		//totalTargetRow("Assigning\n");
		//printf("row:%d col:%d\n",i,i);
		row[i - 27] = i;
		col[i - 27] = i;

	}
//	for (int i = 1; i <= 5; i++) {
//		row[i + 5 - 1] = i + 20;
//		col[i + 5 - 1] = i;
//	}
	for (int i = 0; i < 10; i++) {

		val[i] = 0.001; // (i + 1) * 0.01;
		//printf("val:%lf\n",val[i]);
	}
	for (int i = 0; i < 2; i++) {
		val1[i] = val2[i] = val3[i] = 0.001 * 0.001;
	}
	for (int i = 0; i < 6; i++) {
		hessVal[i] = 0.001 * 0.001;
	}
	for(int i=0;i<2;i++)
	{
		row[i+10]= 39+i;
		col[i+10] = 39+i;
		val[i+10] = 0.01;
	}
//	for(int i=0;i<10;i++)
//	{
//		printf("row:%lf,col:%lf:%lf\n",row[i],col[i],val[i]);
//	}
	influence = new SparseMatrix(row, col, val, 12, 5898240, 1091);
	printf("col:%d\n",influence->getCol());
	printf("Creating hessian\n");
	//hessianM = new SparseMatrix(hessRow, hessCol, hessVal, 6, 1091, 1091, 4);
	hessianM = new SparseMatrix(hessRow, hessCol, hessVal, 6, 1091, 1091);
	hessians[0] = new SparseMatrix(row1, col1, val1, 2, 1091, 1091);
	hessians[1] = new SparseMatrix(row2, col2, val2, 2, 1091, 1091);
	hessians[2] = new SparseMatrix(row3, col3, val3, 2, 1091, 1091);
	//influence = new Matrices(row,col,val,15,500,20);
	printf("col after:%d\n",influence->getCol());
}
/*
 * Get the voxels corresponding to the organ ids provided in the oarId.
 * @param:
 * 	oarId: The organ  id.
 * 	numOrgan: The number of elements in the organ.
 * 	@return the
 */
void oneVoxel::getVoxelofOar(int *oarId,int numOrgan,double *vxls)
{
	int k=0;
	for(int i=0;i<numOrgan;i++)
	{
		for(int j=0;j<numVoxelsRow[oarId[i]];j++)
		{
			vxls[k] = voxels[oarId[i]][j];
			k++;
		}
	}
}
double *oneVoxel::getMetrics(PrioritizedMain::Vector *w1,PrioritizedMain::Vector *w2,PrioritizedMain::Vector *w3,PrioritizedMain::Vector *w4)
{
	double *metrics = new double[getNumTargets()];
	for(int i=0;i<getNumTargets();i++)
	{
		metrics[i] = 75;
	}
	return metrics;
}

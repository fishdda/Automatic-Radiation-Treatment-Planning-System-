/*
 * MatlabCommunicationManager.h
 *
 *  Created on: Nov 7, 2011
 *      Author: user
 */
#include "SparseMatrix.h"
#include "Matrices.h"
#include "engine.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "stdinc.h"
#include "Constant.h"
#ifndef MATLABCOMMUNICATIONMANAGER_H_
#define MATLABCOMMUNICATIONMANAGER_H_

//class MatlabCommunicationManager: public ReferencedObject {
class MatlabCommunicationManager {
	Engine *ep;

	mxArray *targets, *stepIOar, *stepIIOar, *mxStepIIIOar, *bndryVoxels,
			*allVoxels, *mxAlpha, *mxHessRow, *mxHessCol, *mxHessVal,
			*mxRowSize;
	double **dblInfluenceM;

	const mxArray *row, *col, *value;

	int *numVoxels;
	int numTargets;
	int totalTargetRow;
	double ubWeight;
	double *prescDose;
	double maxPrescDose;
	double maxdose;
	int getDimension(int);
	int *stepIorganId;
	int *stepIIorganId;/*Id of the stepII organs.*/
	int *stepIIIorganId; /*Id of the stepII organs.*/
	Matrices *targetMat;
	double *targetsId;
	int ttlStepIOarVoxel;
	int ttlStepIIOarRow;/* Total number of row in the step II organ's voxel.*/
	int numStepIOar;
	int numStepIIOar; /* Number of organs in the step II */
	int numStepIIIOar; /* Number of organs in the step II */
	int ttlStepIIIOarRow;/* Total number of row in the step III organ's voxel.*/
	int nonzeroStepIVox;
	// # of non-zero voxel in step II organ.
	int nonzeroStepIIVoxels;
	int maxOrganId;
	int nonzeroTargetVx;
	double **alpha;
	//double *alpha;
	SparseMatrix **hessians;
	int caseNo;
	void getStepIIVoxels();
	void analyzeAccess();
	Matrices **subMatrix;
	//int nzbmltfortar;
public:
	Matrices **targetInf;
	Matrices **stepIinf;
	SparseMatrix *influence;
	double **voxels;
	SparseMatrix *hessian;
	MatlabCommunicationManager(int caseNum);
	virtual ~MatlabCommunicationManager();
	void breakInfMatrix();
	SparseMatrix *getHMatrix() {
		return hessian;
	}

	SparseMatrix *getHessian(int i) {
		assert(i < getNumTargets());
		return hessians[i];
	}
	int *getStepIOrganId() {
		return stepIorganId;
	}
	int getColInfMat();
	long unsigned int getRowInfMat();

	void getStepIOarVoxels();
	int getNumStepIIOarVoxels();
	int getNumTargetVoxels();
	void analyze(const mxArray *array_ptr);
	void analyze_full(const mxArray *numeric_array_ptr);
	void computerNonZeroVoxelofTargets();
	double* getPrescribedDose();
	double getUBWeight();
	double getMaxDose();
	void getHessians();
	void getStepIIIInfo();
	double getMaxPrescDose();
	void freeDblArray(mxArray *ptr);
	/**
	 * @return the id of organs in the step 2.
	 * @atuthor: Paras Babu Tiwari
	 */
	int * getStepIId() {
		return stepIIorganId;
	}
	/**
	 * @return: the id of organs in step III.
	 * @Author: Paras Babu Tiwari
	 */
	int *getStepIIId() {
		return stepIIIorganId;
	}
	/*
	 * @Input:
	 * 	i: The position of the organ in the array.
	 * @return:
	 * 	The id of the organ.
	 * @Author:
	 * 		Paras Babu Tiwari
	 */
	int getStepIId(int i) {
		return (int) stepIIorganId[i];
	}
	double getDoubleVal(const mxArray *ptr);
	double* getDoubleArray(const mxArray *ptr);
	double ** getDoubleArrays(const mxArray *ptr);
	void getConfigParam();
	void getNzStepIIVoxl();
	void makeSpace();
	//void getTargets();
	void init(int caseNum);
	void getTargetVoxel();
	int getNumTargets();
	void getHessianMatrix();
	double *getMetrics(PrioritizedMain::Vector *w1, PrioritizedMain::Vector *w2,
			PrioritizedMain::Vector *w3,PrioritizedMain::Vector *w4);
	void putCaseId(int caseNum);
	void getOarVoxels(mxArray *oar, double **oarId, int& numOar, int& ttlRow);
	//double *getCellElementAt(const mxArray *cell_ptr,int index);
	//double getNoElement(const mxArray *cell_ptr);

	//double** getAMatrix(double **mat, double *rows, int numRow);

	void getVoxelofOar(int *oarId, int numOrgan, double *vxls);
	double** getVoxels() {
		assert(voxels != NULL);
		return voxels;
	}

	void getInfluenceMatrix();
	/*
	 * @return: Total # of voxels in the target.
	 */
	int getTTLTargetVoxel() {
		//printf("Total target row:%d\n",totalTargetRow);
		assert(totalTargetRow > 0);
		return totalTargetRow;
	}
	int getTargetsId(int i) {
		assert(0 <= i && i <= maxOrganId);
		return (int) targetsId[i];
	}
	/*
	 * @Input:
	 * 	i: The position of the organ in the array.
	 * @return:
	 * 	  The id of the ith organ.
	 * @author: Paras Babu Tiwari
	 */
	int getStepIOrganId(int i) {
		assert(stepIorganId[i] > 0);
		return (int) stepIorganId[i];
	}
	/*
	 * @Input:
	 * 		id: The id of the organ.
	 * @return: The # of voxels of the organ having id i.
	 * @author:Paras Babu Tiwari
	 */
	int* getNumVoxel() {

		return numVoxels;
	}
	/**
	 * @return: Total # of voxels in the stepI organs.
	 */
	int getTTLStepIOar() {
		assert(ttlStepIOarVoxel > 0);
		return ttlStepIOarVoxel;
	}
	int getNumStepIOar() {
		assert(numStepIOar > 0);
		return numStepIOar;
	}
	/*
	 * @return: # of non-zero voxels in the step II organs.
	 * @Author: Paras Babu Tiwari
	 */
	int getNonzeroStepIIVoxels() {
		assert(nonzeroStepIIVoxels > 0);
		return nonzeroStepIIVoxels;

	}
	/*
	 * @return: # of voxels that has non-zero entry in the influence matrix for the step I organ.
	 * The upper bound for this variable is total # of voxels in step I organs * # of beamlets.
	 * @Author: Paras Babu Tiwari
	 */
	int getNonzeroStepIOarVoxels() {

		assert(nonzeroStepIVox > 0);
		return nonzeroStepIVox;
	}
	/*
	 * @return: Total non-zero voxel in the target structure.
	 * @Author: Paras Babu Tiwari
	 */
	int getNonZeroVoxelofTargets() {
		//printf("Total non-zero target voxel:%d\n", nonzeroTargetVx);

		assert(nonzeroTargetVx > 0);
		return nonzeroTargetVx;
	}
	/*Return the number of step II organs.*/
	int getNumStepIIOar() {
		assert(numStepIIOar > 0);
		return numStepIIOar;
	}
	/*
	 * @Input: None
	 * @return the total number of voxels in stepII's organs.
	 * @Author: Paras Babu Tiwari
	 */
	int getTotalStepIIVoxels() {
		int n = 0; //# of voxels
		for (int i = 0; i < getNumStepIIOar(); i++) {
			int id = getStepIId(i);

			n += numVoxels[id];
		}
		assert(n > 0);
		return n;
	}
	/*
	 * @Input: None
	 * @return the total number of voxels in stepII's organs.
	 * @Author: Paras Babu Tiwari
	 */
	int getTotalStepIIIVoxels() {
		int n = 0; //# of voxels
		for (int i = 0; i < numStepIIIOar; i++) {
			int id = stepIIIorganId[i];
			n += numVoxels[id];
		}
		assert(n > 0);
		return n;
	}
	/*Get # of stepIII organs */
	int getNumStepIIIOar() {
		return numStepIIIOar;
	}
	/*
	 * @Param
	 * 	Input:
	 * 		i: Index of the organ in the alpha array.
	 * Return the alpha value
	 */
//	double *getAlpha() {
//
//		return alpha;
//	}
	double **getAlpha() {

		return alpha;
	}
	void getDose(PrioritizedMain::Vector **D, const double *x, int *id,
			int numOar);
	void getDose(PrioritizedMain::Vector **D, PrioritizedMain::Vector *w,
			int *StrId, int numOrgans);
	void getPrescDose();
	int getNumVoxels(int i) {
		return numVoxels[i];
	}
	void getData(double *d,PrioritizedMain::Vector *w);
};

#endif /* MATLABCOMMUNICATIONMANAGER_H_ */

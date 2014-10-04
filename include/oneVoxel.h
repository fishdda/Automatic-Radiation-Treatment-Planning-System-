/*
 * oneVoxel.h
 *
 *  Created on: Dec 28, 2011
 *      Author: user
 */
#include "FileManager.h"
#include "SparseMatrix.h"
#include "stdinc.h"

#ifndef ONEVOXEL_H_
#define ONEVOXEL_H_

class oneVoxel: public ReferencedObject {
public:
	oneVoxel();
	virtual ~oneVoxel();
private:
	int* numVoxelsRow;
	int numTargets;
	int *stepIorganId;
	int *stepIIorganId;
	int *stepIIIorganId;
	int *targetsId;
	int numStepIOar;
	int numStepIIIOar;
	int ttlTargetRow;
	int ttlStepIOarRow;
	double *prescDose;
	double **alpha;
	SparseMatrix *hessianM;
	SparseMatrix **hessians;
public:
	SparseMatrix *getHMatrix() {
		return hessianM;
	}
	SparseMatrix *influence;
	double **voxels;
	SparseMatrix *getHessian(int i) {
		assert(i < getNumTargets());
		return hessians[i];
	}
	int getNonzeroStepIIVoxels() {
		return 2;
	}
	void getStepIIVoxel();
	void getStepIIIinfo();
	int* getNumVoxel()
	{
		return numVoxelsRow;
	}
	int getNumVoxels(int i) {
		return numVoxelsRow[i];
	}
	int getNumTargets() {
		return numTargets;
	}
	int *getStepIOrganId() {
		return stepIorganId;
	}
	int getStepIOrganId(int i) {
			return stepIorganId[i];
		}
	int getTargetsId(int i) {
		return targetsId[i];
	}
	int *getTargetsId()
	{
		return targetsId;
	}
	int getNumStepIOar() {
		return numStepIOar;
	}
	int getTTLTargetVoxel() {
		return ttlTargetRow;
	}
	int getTTLStepIOar() {
		return ttlStepIOarRow;
	}

	void getTargetVoxel();
	void getInfluenceMatrix();
	void init();
	int getColInfMat() {

//		printf("influence:%x\n",influence);
//		printf("%d\n",influence->getCol());
		if(influence != NULL)
			return influence->getCol();
		else return 0;
	}
	long unsigned int getRowInfMat() {
		return influence->getRow();
	}

	int getnumTargets() {
		return numTargets;
	}
	double** getVoxels() {
		return voxels;
	}

	double* getPrescribedDose() {

		return prescDose;
	}
	double getMaxDose() {
		return 80;
	}
	double getUBWeight() {
		return 42.7559;
	}

	int getNonzeroStepIOarVoxels() {

		return 2;
	}

	int getNonZeroVoxelofTargets() {
		return 6;
		///return 1;
	}
	/*Return the number of step II organs.*/
	int getNumStepIIOar() {
		return 1;
	}
	int getTotalStepIIVoxels() {
		return numVoxelsRow[stepIIorganId[0]];
	}
	int getNumStepIIIOar() {
		return numStepIIIOar;
	}
	/*
	 * @Param
	 * 	Input:
	 * 		i: Index of the organ in the alpha array.
	 * Return the alpha value
	 */
	double **getAlpha() {

		return alpha;
	}
	void getVoxelStepI();
	void getDose(PrioritizedMain::Vector **D, const double *x, int *id, int numOar);
	void getDose(PrioritizedMain::Vector **D,PrioritizedMain::Vector *w, int *StrId, int numOrgans);
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
	int * getStepIId() {
		return  stepIIorganId;
	}
	int *getStepIIId() {
		return stepIIIorganId;
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
			n += numVoxelsRow[id];
			//printf("id:%d\n",id);
		}
		assert(n > 0);
		return n;
	}
	void getVoxelofOar(int *oarId,int numOrgan,double *vxls);
	double *getMetrics(PrioritizedMain::Vector *w1,PrioritizedMain::Vector *w2,PrioritizedMain::Vector *w3,PrioritizedMain::Vector *w4);

};

#endif /* ONEVOXEL_H_ */

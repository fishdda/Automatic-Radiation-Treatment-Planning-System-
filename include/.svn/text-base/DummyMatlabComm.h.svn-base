/*
 * DummyMatlabComm.h
 *
 *  Created on: Nov 23, 2011
 *      Author: user
 */
#include <stdio.h>
//#include "Matrix.h"
#include "SparseMatrix.h"
#ifndef DUMMYMATLABCOMM_H_
#define DUMMYMATLABCOMM_H_

class DummyMatlabComm {
private:
	int *numVoxelsRow;
	int numTargets;
	int *stepIorganId;
	int *targetsId;
	int numStepIOar;
	int ttlTargetRow;
	int ttlStepIOarRow;
protected:
	SparseMatrix *influence;
	double **voxels;
public:
	int getNumVoxelsRow(int i) {
		return numVoxelsRow[i];
	}
	int getNumTargets() {
		return numTargets;
	}
	int getStepIOrganId(int i) {
		return stepIorganId[i];
	}
	int getTargetsId(int i) {
		return targetsId[i];
	}
	int getNumStepIOar() {
		return numStepIOar;
	}
	int getTTLTargetRow()
	{
			return ttlTargetRow;
	}
	int getTTLStepIOar()
	{
		return ttlStepIOarRow;
	}
	DummyMatlabComm();
	virtual ~DummyMatlabComm();

	void getTargetVoxel();
	void getInfluenceMatrix();
	void init();
	int getColInfMat() {

		return influence->getCol();
	}
	int getRowInfMat() {
		return influence->getRow();
	}

	int getnumTargets() {
		return numTargets;
	}
	double** getVoxels() {
		return voxels;
	}

	double getPrescribedDose() {

		return 70;
	}
	double getMaxDose() {
		return 80;
	}
	double getUBWeight() {
		return 42.7559;
	}
	int getNumStepI() {

		return 1;
	}
	int getNonzeroStepIOarVoxels() {

		return 105142;
	}

	int getNonZeroVoxelofTargets() {
		return 2816997;
		//return getTTLTargetRow() * getColInfMat();
	}
	void getVoxelStepI();
};

#endif /* DUMMYMATLABCOMM_H_ */

/*
 * Result.h
 *
 *  Created on: Jul 17, 2012
 *      Author: user
 */

#ifndef RESULT_H_
#define RESULT_H_
class StepImpl;
class StepIImpl;
#include <assert.h>
#include "Matrices.h"

class Result {
	//Beamlet Found in the Step I.
	Matrices *w1;
	//Beamlet Found in the Step II.
	Matrices *w2;
	//Pointer to the StepI implementation
	StepImpl *stepI;
	StepIImpl *stepII;
	//Beamlet Size
	int bmltSize;

	//Current Step #;
	int step;

public:
	Result();
	virtual ~Result();

	void setw1(Matrices *w) {

		this->w1 = w;
	}
	void setw2(Matrices *w) {
		this->w2 = w;
	}
	void setBlmtSize(int s) {
		bmltSize = s;
	}
	void setStepI(StepImpl *s) {
		stepI = s;
	}
	void PrintDose();
	void setStepPtr(void *ptr);
	void setStep(int s){step = s;}
	void setW(const double *x);
	Matrices *getW1(){return w1;}
	Matrices *getW2(){return w2;}
	StepImpl *StepI(){return stepI;}

};

#endif /* RESULT_H_ */

/*
 * Result.cpp
 *
 *  Created on: Jul 17, 2012
 *      Author: user
 */

#include "Result.h"

Result::Result() {
	// TODO Auto-generated constructor stub
	stepI = NULL;

}

Result::~Result() {
	// TODO Auto-generated destructor stub
}

void Result::setStepPtr(void *stpPtr) {
	assert(stpPtr != NULL);
	printf("step:%d\n",step);
	switch (step) {
		case 1:
			stepI = (StepImpl *) stpPtr;
			break;
		case 2:
			stepII = (StepIImpl *) stpPtr;
			break;
	}
}
void Result::setW(const double *x) {
	assert(bmltSize > 0);
	Matrices *w;
	//Matrices *w = new Matrices(x, bmltSize);
	//printf("w\n");
	//w->printNonZeroEntry();
	switch (step) {
		case 1:
			setw1(w);
			break;
		case 2:
			this->w2 = w;
			break;
	}
}

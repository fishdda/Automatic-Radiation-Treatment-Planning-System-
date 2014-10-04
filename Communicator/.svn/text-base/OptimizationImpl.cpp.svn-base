/*
 * OptimizationImpl.cpp
 *
 *  Created on: Jul 11, 2012
 *      Author: user
 */

#include "OptimizationImpl.h"

OptimizationImpl::OptimizationImpl(int step, ProblemReps *prob) {
	// TODO Auto-generated constructor stub
	//Call the MatlabCommunication Manager init function
	this->step = step;
	current = prob;

}

OptimizationImpl::~OptimizationImpl() {
	// TODO Auto-generated destructor stub

}
/*
 * Set the current step of the optimization;
 * @Param
 * 	Input
 * 		s: Step Number
 * 	@author: Paras Babu Tiwari
 *
 */
//void OptimizationImpl::setStep(int s, ProblemReps *prob) {
//
//}
/**
 * Evaluate the objective function
 */
long double OptimizationImpl::evaluateObjective(const double *x) {
	return current->evaluateObjective(x);

}
void OptimizationImpl::getStartingPoint(double *x,double *z_l,double *z_u) {
	return current->getStartingPoint(x,z_l,z_u);

}
/*
 * Evaluate the gradient of the function.
 */
void OptimizationImpl::evaluategradf(const double *x, double *g) {
	current->evaluategradf(x, g);

}
/**
 * Evaulate the constraints
 */
void OptimizationImpl::evalg(const double *x, double *g) {

	assert(current != NULL);

	current->evalg(x, g);

}
/**
 * Get the number of variables.
 */
int OptimizationImpl::getNumVar() {
	//printf("stepI:%x\n",stepI);
	//printf("Step:%d\n",step);

	current->getNumVar();

}
/**
 * Get the number of constraints.
 */
int OptimizationImpl::getNumConstraint() {

	current->getNumConstraint();

}
/**
 * Get # of non-zero Jacobian entries.
 */
int OptimizationImpl::getNumJac() {
	current->getNumJac();

}
/**
 * Get the Hessaian Matrices.
 */
bool OptimizationImpl::getHessian(int *iRow, int *iCol, double *value,
		double obj_factor,const double *lambda ) {
	return current->getHessian(iRow,iCol,value,obj_factor,lambda);

}
/**
 * Get the number of non-zero hessian entry
 */
int OptimizationImpl::getNumHess() {
	current->getNumHess();

}
/**
 * Get the Jacobian Matrices
 */
void OptimizationImpl::getJacobian(int * iRow, int * iCol, double* value, const double *x) {

	current->getJacobian(iRow,iCol,value,x);

}
/**
 * Get the Bounds info.
 * Input Param: None
 * @Output Param:
 * 	x_l: Lower bound of variables.
 * 	x_u: Upper bound of variables.
 */
void OptimizationImpl::getBounds(double *x_l, double *x_u, double *g_l,
		double *g_u) {

	current->getBounds(x_l,x_u,g_l,g_u);

}
//void OptimizationImpl::printDose(Matrices *w) {
//
//	current->printDose(w);
//
//}
///**
// * Set the w matrix obtained by the optimization.
// * @param w: 'w' matrix obtained from the optimization.
// * @return : None
// */
//int OptimizationImpl::setw(Matrices *w) {
//
//	current->setw(x, g);
//
//}
/**
 * Print Dose
 */
void OptimizationImpl::PrintDose() {

	current->printDose();

}
/*
 * Set the Objective Function Value
 */
void OptimizationImpl::setObjectiveFunValue(double val)
{
	current->setObjectiveFunValue(val);
}
/**
 *
 */
void OptimizationImpl::setW(const double *x,const double *z_l,const double *z_u,const double *lambda)
{
	current->setW(x,z_l,z_u,lambda);
}

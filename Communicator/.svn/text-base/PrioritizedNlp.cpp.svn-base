// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: PrioritizedNlp.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-16

#include "PrioritizedNlp.hpp"

using namespace std;
using namespace Ipopt;

// constructor
PrioritizedNlp::PrioritizedNlp(int step, ProblemReps *prob) {
	//opt = new MatlabCommunicationManager();
	//opt = new DummyMatlabComm();
	opt = new OptimizationImpl(step, prob);
	fdiff=0;
	grdiff=0;
	gdiff=0;
	hdiff=0;
	jacdiff=0;
}

//destructor
PrioritizedNlp::~PrioritizedNlp() {
}
/*
 * Set the step # and the pointer to the current step.
 * @param:
 * 	Input:
 * 		step: Current step #.
 * 		prob: Pointer to the object of current step.
 * @author: Paras Babu Tiwari
 */
//void PrioritizedNlp::init(int step,ProblemReps *prob)
//{
//	opt->setStep(step,prob);
//}
// returns the size of the problem
bool PrioritizedNlp::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
		Index& nnz_h_lag, IndexStyleEnum& index_style) {

	//printf("Inside problem information\n");
	n = opt->getNumVar();
	printf("# of variable:%d\n", n);
	m = opt->getNumConstraint();
	printf("# of constraints:%d\n", m);
	nnz_jac_g = opt->getNumJac();
	printf("# of Jacobian:%d\n", nnz_jac_g);
	//nnz_h_lag = 0;
	nnz_h_lag = opt->getNumHess();
	printf("# of Hessian:%d\n", nnz_h_lag);
	// use the C style indexing (0-based)
	index_style = TNLP::C_STYLE;
	//printf("Done getting problem information\n");
	return true;
}

// returns the variable bounds
bool PrioritizedNlp::get_bounds_info(Index n, Number* x_l, Number* x_u,
		Index m, Number* g_l, Number* g_u) {
	// here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
	// If desired, we could assert to make sure they are what we think they are.
	int totalvar = opt->getNumVar();
	int totalconst = opt->getNumConstraint();
	assert(n == totalvar);
	assert(m == totalconst);
	int i = 0;
	//printf("Inside get_bounds_info\n");
	opt->getBounds(x_l, x_u, g_l, g_u);
	//printf("Done Getting bound info\n");
	return true;
}

// returns the initial point for the problem
bool PrioritizedNlp::get_starting_point(Index n, bool init_x, Number* x,
		bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda,
		Number* lambda) {
	// Here, we assume we only have starting values for x, if you code
	// your own NLP, you can provide starting values for the dual variables
	// if you wish
	assert(init_x == true);
	assert(init_z == false);
	assert(init_lambda == false);
	//printf("Getting starting point \n");
	//printf("Inside Starting point\n");
	// initialize to the given starting point
	//for (Index i = 0; i < opt->getNumVar(); i++) {
	//	x[i] = 2;
	//}
	opt->getStartingPoint(x,z_L,z_U);
	//printf("Done Starting point\n");
	return true;
}

// returns the value of the objective function
bool PrioritizedNlp::eval_f(Index n, const Number* x, bool new_x,
		Number& obj_value) {
	//printf("Inside eval f\n");
	int totalVar = 0;
	totalVar = opt->getNumVar();
	assert(n == totalVar);
	//printf("Inside Objective function\n");
	time_t start_t = time(NULL);
	obj_value = opt->evaluateObjective(x);
	time_t end_t = time(NULL);
	double diff = difftime(end_t, start_t);
	fdiff = fdiff + diff;
	//printf("Done Objective function\n");

	//printf("Objective function:%lf\n", obj_value);
	return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool PrioritizedNlp::eval_grad_f(Index n, const Number* x, bool new_x,
		Number* grad_f) {
	//printf("Inside gradf");
	time_t start_t = time(NULL);
	int totalVar = opt->getNumVar();
	assert(n == totalVar);
	//printf("Evaluate grad\n");
	double *g;
	//printf("Inside grad_f\n");
	int s = opt->getNumVar();
	g = new double[s];
	opt->evaluategradf(x, g);

	for (int i = 0; i < totalVar; i++) {
		//		printf("grad_f[%d]:%Lf\n",i,g[i]);
		grad_f[i] = g[i];

	}
	time_t end_t = time(NULL);
	double diff = difftime(end_t, start_t);
	grdiff = grdiff + diff;
	//printf("Gradient time:%lf\n",grdiff);
	delete[] g;
	g = NULL;
	//printf("Done grad_f\n");
	return true;
}

// return the value of the constraints: g(x)
bool PrioritizedNlp::eval_g(Index n, const Number* x, bool new_x, Index m,
		Number* g) {
	time_t start_t = time(NULL);
	//printf("Eval_g\n");
	assert(n == opt->getNumVar());
	assert(m == opt->getNumConstraint());
	//printf("Inside eval_g\n");

	opt->evalg(x, g);
	//printf("Done eval_g\n");
	time_t end_t = time(NULL);
	double diff = difftime(end_t, start_t);
	gdiff = gdiff + diff;

	return true;
}

// return the structure or values of the jacobian
bool PrioritizedNlp::eval_jac_g(Index n, const Number* x, bool new_x, Index m,
		Index nele_jac, Index* iRow, Index *jCol, Number* values) {
	//printf("Getting Jacobian\n");
	//printf("Inside Jacobian\n");
	time_t start_t = time(NULL);
	opt->getJacobian(iRow, jCol, values, x);
	time_t end_t = time(NULL);
	double diff = difftime(end_t, start_t);
	jacdiff = jacdiff + diff;
	//printf("Done With Jacobian\n");
	return true;
}

//return the structure or values of the hessian
bool PrioritizedNlp::eval_h(Index n, const Number* x, bool new_x,
		Number obj_factor, Index m, const Number* lambda, bool new_lambda,
		Index nele_hess, Index* iRow, Index* jCol, Number* values) {

	//printf("Inside hessian\n");
	time_t start_t = time(NULL);
	return opt->getHessian(iRow, jCol, values, obj_factor, lambda);
	time_t end_t = time(NULL);
	double diff = difftime(end_t, start_t);
	hdiff = hdiff + diff;
	//printf("Done with hessian\n");
}
void PrioritizedNlp::finalize_solution(SolverReturn status, Index n,
		const Number* x, const Number* z_L, const Number* z_U, Index m,
		const Number* g, const Number* lambda, Number obj_value,
		const IpoptData* ip_data, IpoptCalculatedQuantities* ip_cq) {

	//#ifdef dummy
	//	 t = ip->TimingStats().OverallAlgorithm().TotalCpuTime();
	//#else
	//	 t = ip->TimingStats().OverallAlgorithm().TotalTime();
	//
	//#endif
	//	t = ip->TimingStats().OverallAlgorithm().TotalCpuTime();
	//	//cout<<"CPU Time:"<<t<<endl;
	//	printf("Time:%lf\n",t);
	//	f->writeFile(t);
	//double t = ip_data->TimingStats().OverallAlgorithm().TotalCpuTime();
	printf("Time in Objective:%lf\n", fdiff);
	printf("Time in Gradient:%lf\n", grdiff);
	printf("Time in Constraint:%lf\n", gdiff);
	printf("Time in Hessian:%lf\n", hdiff);
	printf("Time in Jacobian:%lf\n", jacdiff);
	printf("\n\nSolution of the primal variables, x\n");
	opt->setW(x,z_L,z_U,lambda);
	opt->PrintDose();
	opt->setObjectiveFunValue(obj_value);
}
void PrioritizedNlp::finalize_solution(SolverReturn status, Index n,
		const Number* x, const Number* z_L, const Number* z_U, Index m,
		const Number* g, const Number* lambda, Number obj_value) {
	// here is where we would store the solution to variables, or write to a file, etc
	// so we could use the solution.

	opt->setW(x,z_L,z_U,lambda);
	opt->PrintDose();
	opt->setObjectiveFunValue(obj_value);
}

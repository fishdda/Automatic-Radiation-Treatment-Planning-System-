// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: hs071_main.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-10

#include "StepIImpl.h"
#include "StepImpl.h"
#include "StepIIImpl.h"
#include "StepIVmpl.h"
#include "stdinc.h"
#ifdef dummy
#include "oneVoxel.h"
#else
#include "MatlabCommunicationManager.h"
#endif
#include "PrioritizedNlp.hpp"
//#include "IpIpoptApplication.hpp"
#include "Constant.h"
//using namespace Ipopt;
#define NOSTPS 4

int main(int argv, char* argc[]) {
	printf("inside the main\n");
	char buffer[50];
	if(argv<2){
		printf("Invalid command. Should be ./PrioritizedNlp <fileName>\n");
		exit(0);
	}
	printf("Reading content from:%s\n",argc[1]);
	FileManager fp(argc[1], Constant::Read);

	//printf("Constant:%lf\n",Constant::NUMCASES);
	double t =0,total=0;
	std::vector<int> cases;
	fp.readFile(cases);
	fp.close();
	for (int i = 0; i < cases.size(); i++) {
		//MatlabCommunicationManager *mtlb = new MatlabCommunicationManager(Constant::CASEID[i]);
#ifdef dummy
		//SmartPtr<oneVoxel> mtlb = new oneVoxel();
		oneVoxel *mtlb = new oneVoxel();
#else
//		SmartPtr<MatlabCommunicationManager> mtlb =
//				new MatlabCommunicationManager(Constant::CASEID[i]);
		printf("Case:%d\n",cases[i]);
		MatlabCommunicationManager *mtlb = new MatlabCommunicationManager(
				cases[i]);
#endif
		stringstream ss;
		ss << cases[i];
		string str = ss.str();
		std::string fname = "time" + str+".txt";
		const char* nm = fname.c_str();
		FileManager ftime(nm, Constant::Write);
		//printf("col:%d\n",mtlb->getColInfMat());
		//printf("Case Number:%d\n",Constant::CASEID[i]);
		// Create a new instance of your nlp
		//  (use a SmartPtr, not raw)

		// Create a new instance of IpoptApplication
		//  (use a SmartPtr, not raw)
		//SmartPtr < IpoptApplication > app = new IpoptApplication();

		// Change some options
		// Note: The following choices are only examples, they might not be
		//       suitable for your optimization problem.
		time_t start_t, end_t;
		start_t = time(NULL);
		printf("%ld hours since January 1, 1970\n", start_t);

		// Intialize the IpoptApplication and process the options
		SmartPtr<TNLP> mynlp;
		ApplicationReturnStatus status;
		StepImpl *s1;
		StepIImpl *s2;
		StepIIImpl *s3;
		StepIVmpl *s4;
		//oneVoxel *mtlb = new oneVoxel();
		for (int j = 0; j < NOSTPS; j++) {
			if (j == 0) {

				SmartPtr<IpoptApplication> app = new IpoptApplication();
				app->Options()->SetStringValue("mu_strategy", "adaptive");
				app->Options()->SetStringValue("linear_solver","ma86");
				//app->Options()->SetStringValue("output_file", "ipopt.out");
				//app->Options()->SetStringValue("derivative_test","second-order");
				//app->Options()->SetStringValue("hessian_approximation", "limited-memory");
				app->Initialize();

				s1 = new StepImpl(mtlb);
				mynlp = new PrioritizedNlp(i + 1, s1);
				printf("Starting Optimization\n");
				status = app->OptimizeTNLP(mynlp);
				t = app->IpoptDataObject()->TimingStats().OverallAlgorithm().TotalCpuTime();
				total+=t;
				//ftime.writeFile(t);

			} else if (j == 1) {
				SmartPtr<IpoptApplication> app = new IpoptApplication();
				app->Options()->SetStringValue("mu_strategy", "adaptive");
				app->Options()->SetStringValue("output_file", "ipopt1.out");
				app->Options()->SetStringValue("linear_solver","ma86");
				//app->Options()->SetNumericValue("warm_start_bound_push",1e-6);
				//app->Options()->SetStringValue("warm_start_init_point","yes");
				//app->Options()->SetStringValue("derivative_test",	"second-order");
				app->Initialize();
				s2 = new StepIImpl(mtlb, s1);
				mynlp = new PrioritizedNlp(i + 1, s2);
				status = app->OptimizeTNLP(mynlp);
				t = app->IpoptDataObject()->TimingStats().OverallAlgorithm().TotalCpuTime();
				total+=t;
				//ftime.writeFile(t);
			} else if (j == 2) {

				SmartPtr<IpoptApplication> app = new IpoptApplication();
				app->Options()->SetStringValue("mu_strategy", "adaptive");
				app->Options()->SetStringValue("linear_solver","ma86");
				//app->Options()->SetStringValue("output_file", buffer);
				//app->Options()->SetStringValue("hessian_approximation", "limited-memory");
				//app->Options()->SetStringValue("derivative_test","second-order");
				app->Initialize();
				printf("just before col:%d\n", mtlb->getColInfMat());
				s3 = new StepIIImpl(mtlb, s2);
				mynlp = new PrioritizedNlp(j + 1, s3);
				status = app->OptimizeTNLP(mynlp);
				t = app->IpoptDataObject()->TimingStats().OverallAlgorithm().TotalCpuTime();
				total+=t;
				//ftime.writeFile(t);
			} else if (j == 3) {

				SmartPtr<IpoptApplication> app = new IpoptApplication();
				app->Options()->SetStringValue("mu_strategy", "adaptive");
				app->Options()->SetStringValue("linear_solver","ma86");
				//app->Options()->SetStringValue("output_file", buffer);
				//app->Options()->SetStringValue("hessian_approximation", "limited-memory");
				//app->Options()->SetStringValue("derivative_test","second-order");
				app->Initialize();
				//printf("just before col:%d\n", mtlb->getColInfMat());
				s4 = new StepIVmpl(mtlb, s3);
				mynlp = new PrioritizedNlp(j + 1, s4);
				status = app->OptimizeTNLP(mynlp);
				t = app->IpoptDataObject()->TimingStats().OverallAlgorithm().TotalCpuTime();
				total+=t;
				//ftime.writeFile(t);
			}
			//mynlp->init(i+1,mtlb);
			// Ask Ipopt to solve the problem
			//ApplicationReturnStatus status = app->OptimizeTNLP(mynlp);
		}
		ftime.writeFile(total);

		end_t = time(NULL);
		printf("%ld hours since January 1, 1970", end_t);
		double diff = difftime(end_t, start_t);
		printf("Elapsed time:%lf seconds\n", diff);
		//delete mtlb;
		//cout << "Time elapsed: " << double(diffclock(end,begin)) << " ms"<< endl;
		// As the SmartPtrs go out of scope, the reference count
		// will be decremented and the objects will automatically
		// be deleted.
		ftime.close();
	}

	return (int) 0;
}


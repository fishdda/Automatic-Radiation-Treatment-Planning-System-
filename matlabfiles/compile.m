cd .

mex  -f engopts.sh -I/Users/parastiwari/Ipopt3103/include/coin -I/Users/parastiwari/Documents/workspace/Ipopt/src/include -I/usr/local/include/ -c ../Communicator/PrioritizedNlp.cpp;
mex  -f engopts.sh -I/Users/parastiwari/Ipopt3103/include/coin -I/Users/parastiwari/Documents/workspace/Ipopt/src/include -I/usr/local/include/ -c ../Communicator/PrioritizedMain.cpp;
mex  -f engopts.sh -I/Users/parastiwari/Ipopt3103/include/coin -I/Users/parastiwari/Documents/workspace/Ipopt/src/include -I/usr/local/include/ -c ../Communicator/MatlabCommunicationManager.cpp;
mex  -f engopts.sh -I/Users/parastiwari/Ipopt3103/include/coin -I/Users/parastiwari/Documents/workspace/Ipopt/src/include -I/usr/local/include/ -c ../Communicator/OptimizationImpl.cpp;
mex  -f engopts.sh -I/Users/parastiwari/Ipopt3103/include/coin -I/Users/parastiwari/Documents/workspace/Ipopt/src/include -I/usr/local/include/ -c ../datastructure/Matrices.cpp;
mex  -f engopts.sh -I/Users/parastiwari/Ipopt3103/include/coin -I/Users/parastiwari/Documents/workspace/Ipopt/src/include -I/usr/local/include/ -c ../datastructure/SparseMatrix.cpp;
mex  -f engopts.sh -I/Users/parastiwari/Ipopt3103/include/coin -I/Users/parastiwari/Documents/workspace/Ipopt/src/include -I/usr/local/include/ -c ../datastructure/Vector.cpp;
mex  -f engopts.sh -I/Users/parastiwari/Ipopt3103/include/coin -I/Users/parastiwari/Documents/workspace/Ipopt/src/include -I/usr/local/include/ -c ../Communicator/Util.cpp;
mex  -f engopts.sh -I/Users/parastiwari/Ipopt3103/include/coin -I/Users/parastiwari/Documents/workspace/Ipopt/src/include -I/usr/local/include/ -c ../Communicator/FileManager.cpp;
mex  -f engopts.sh -I/Users/parastiwari/Ipopt3103/include/coin -I/Users/parastiwari/Documents/workspace/Ipopt/src/include -I/usr/local/include/ -c ../steps/StepImpl.cpp;
mex  -f engopts.sh -I/Users/parastiwari/Ipopt3103/include/coin -I/Users/parastiwari/Documents/workspace/Ipopt/src/include -I/usr/local/include/ -c ../steps/StepIImpl.cpp;
mex  -f engopts.sh -I/Users/parastiwari/Ipopt3103/include/coin -I/Users/parastiwari/Documents/workspace/Ipopt/src/include -I/usr/local/include/ -c ../steps/StepIIImpl.cpp;
mex  -f engopts.sh -I/Users/parastiwari/Ipopt3103/include/coin -I/Users/parastiwari/Documents/workspace/Ipopt/src/include -I/usr/local/include/ -c ../steps/StepIVmpl.cpp;
mex  -f engopts.sh -I/Users/parastiwari/Ipopt3103/include/coin -I/Users/parastiwari/Documents/workspace/Ipopt/src/include -I/usr/local/include/ -c ../Communicator/ProblemReps.cpp;

mex  -f engopts.sh  -I/usr/local/include/   -I/Users/parastiwari/Ipopt3103/include/coin -I/Users/parastiwari/Documents/workspace/Ipopt/src/include ../Communicator/PrioritizedNlp.cpp ../steps/StepIImpl.cpp ../steps/StepIIImpl.cpp ../Communicator/PrioritizedMain.cpp ../Communicator/MatlabCommunicationManager.cpp ../Communicator/OptimizationImpl.cpp ../steps/StepIVmpl.cpp ../datastructure/Matrices.cpp ../datastructure/SparseMatrix.cpp ../datastructure/Vector.cpp  ../Communicator/Util.cpp ../Communicator/FileManager.cpp ../steps/StepImpl.cpp  ../Communicator/ProblemReps.cpp  -L/Users/parastiwari/Ipopt3103/lib -lcoinhsl -lipopt -llapack -lblas -lm -lcoinblas;



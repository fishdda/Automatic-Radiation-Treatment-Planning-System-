cd .

mex  -f engopts.sh -I/home/research1/tiwarip/CoinIpopt/include/ipopt -c PrioritizedNlp.cpp
mex  -f engopts.sh -I/home/research1/tiwarip/CoinIpopt/include/ipopt -c PrioritizedMain.cpp
mex  -f engopts.sh -I/home/research1/tiwarip/CoinIpopt/include/ipopt -c oneVoxel.cpp
mex  -f engopts.sh -I/home/research1/tiwarip/CoinIpopt/include/ipopt -c OptimizationImpl.cpp
mex  -f engopts.sh -I/home/research1/tiwarip/CoinIpopt/include/ipopt -c Matrix.cpp
mex  -f engopts.sh  I-/home/research1/tiwarip//CoinIpopt/include/ipopt -c SparseMatrix.cpp
mex  -f engopts.sh -I/home/research1/tiwarip/CoinIpopt/include/ipopt -c HashTable.cpp
mex  -f engopts.sh -I/home/research1/tiwarip/CoinIpopt/include/ipopt -c HashMap.cpp
mex  -f engopts.sh -I/home/research1/tiwarip/CoinIpopt/include/ipopt -c Util.cpp
mex  -f engopts.sh -I/home/research1/tiwarip/CoinIpopt/include/ipopt -c FileManager.cpp
mex  -f engopts.sh -I/home/research1/tiwarip/CoinIpopt/include/ipopt -c StepIIImpl.cpp
mex  -f engopts.sh -I/home/research1/tiwarip/CoinIpopt/include/ipopt -c ProblemReps.cpp
mex  -f engopts.sh -I/home/research1/tiwarip/CoinIpopt/include/ipopt -c UiSetPair.cpp
mex  -f engopts.sh -I/home/research1/tiwarip/CoinIpopt/include/ipopt -c Constant.cpp

mex  -f engopts.sh  -g  -I/home/research1/tiwarip/CoinIpopt/include/ipopt  PrioritizedMain.cpp Constant.cpp UiSetPair.cpp StepIIImpl.cpp ProblemReps.cpp Matrix.cpp ProblemReps.h Util.cpp FileManager.cpp HashMap.cpp SparseMatrix.cpp HashTable.cpp   PrioritizedNlp.cpp oneVoxel.cpp OptimizationImpl.cpp  -L/home/research1/tiwarip/CoinIpopt/lib -lipopt -llapack -lblas -lm



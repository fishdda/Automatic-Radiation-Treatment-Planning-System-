cd .
mex  -f engopts.sh -I/home/research1/tiwarip/CoinIpopt/include/ipopt -c hs071_nlp.cpp
mex  -f engopts.sh -I/home/research1/tiwarip/CoinIpopt/include/ipopt -c hs071_main.cpp
mex  -f engopts.sh -I/home/research1/tiwarip/CoinIpopt/include/ipopt -c DummyMatlabComm.cpp
mex  -f engopts.sh -I/home/research1/tiwarip/CoinIpopt/include/ipopt -c OptFormImpl.cpp
mex  -f engopts.sh -I/home/research1/tiwarip/CoinIpopt/include/ipopt -c Matrix.cpp
mex  -f engopts.sh -I/home/research1/tiwarip/CoinIpopt/include/ipopt -c SparseMatrix.cpp
mex  -f engopts.sh -I/home/research1/tiwarip/CoinIpopt/include/ipopt -c HashTable.cpp
mex  -f engopts.sh -I/home/research1/tiwarip/CoinIpopt/include/ipopt -c FileManager.cpp

mex  -f engopts.sh   -I/home/research1/tiwarip/CoinIpopt/include/ipopt   hs071_main.cpp Util.cpp HashMap.cpp  FileManager.cpp SparseMatrix.cpp HashTable.cpp Matrix.cpp  hs071_nlp.cpp DummyMatlabComm.cpp OptFormImpl.cpp  -L/home/research1/tiwarip/CoinIpopt/lib -lipopt -llapack -lblas -lm



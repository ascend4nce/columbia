BNAME?=default
TNAME?=simple

CPP = icpc
OPTS = -O3
MKL = -mkl=sequential

builder:
	$(CPP) builder.cpp newton.cpp maxvol.cpp columbia.cpp nonlinearities.cpp fullscheme.cpp matrixCSR.cpp $(MKL) $(OPTS) -o force71
	#sbatch runbuilder.sh $(BNAME)

simple:
	$(CPP) simpletest.cpp newton.cpp columbia.cpp nonlinearities.cpp fullscheme.cpp matrixCSR.cpp maxvol.cpp $(MKL) $(OPTS) -o force72
	#sbatch runsimple.sh $(TNAME)

reduceone:
	$(CPP) reduced.cpp newton.cpp maxvol.cpp columbia.cpp nonlinearities.cpp fullscheme.cpp reduction.cpp matrixCSR.cpp $(MKL) $(OPTS) -o force73
	#sbatch runreduceone.sh $(BNAME)

reducefull:
	$(CPP) reducedfull.cpp newton.cpp maxvol.cpp columbia.cpp nonlinearities.cpp fullscheme.cpp reduction.cpp matrixCSR.cpp $(MKL) $(OPTS) -o force76
	#sbatch runreducefull.sh $(BNAME)

baseprint:
	$(CPP) printbasis.cpp newton.cpp maxvol.cpp columbia.cpp nonlinearities.cpp fullscheme.cpp reduction.cpp matrixCSR.cpp $(MKL) $(OPTS) -o force75
	#sbatch runbaseprint.sh $(BNAME)

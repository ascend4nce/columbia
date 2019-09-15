#include <mkl.h>
#include <mkl_rci.h>
#include <mkl_blas.h>
#include <mkl_spblas.h>
#include <mkl_service.h>
#include <cmath>
#include <cstdio>
#include "matrixCSR.h"
#include "nonlinearities.h"
#include "fullscheme.h"
#include "reduction.h"
#include "columbia.h"
#include "newton.h"
#include "constants.h"

int main(int argc, char ** argv)
{	
	double * point, * refpoint;
	int gridN;
	double shift, alpha, beta;

	char namev[60];
	char namefx[60];
	char namefy[60];
	char nameB[60];
	char nameH[60];

	if (argc < 2)
	{
		return -1;
	};

	solverparams PRS(GRIDsize, Alpha, Beta, Shiftx, Shifty, Scaling);

	point = new double[PRS.varnum];
	refpoint = new double[PRS.varnum];

	sprintf(namev, "stest/%svalue.dat", argv[1]);
	sprintf(namefx, "stest/%sfluxx.dat", argv[1]);
	sprintf(namefy, "stest/%sfluxy.dat", argv[1]);
	sprintf(nameB, "stest/%smoduleB.dat", argv[1]);
	sprintf(nameH, "stest/%smoduleH.dat", argv[1]);
	
	//simpletest(&PRS, point);
	multitest(&PRS, point, GridNumber);	

	//putsolution(&PRS, refpoint);
	//printrs(&PRS, "stest/rs.dat");

	printsimple(&PRS, point, namev, namefx, namefy);	
	printfluxmodule(&PRS, point, nameH, nameB);
	//printreference(&PRS, point, refpoint, "stest/exp");

	delete[] point;
	delete[] refpoint;
};

#include <mkl.h>
#include <mkl_rci.h>
#include <mkl_blas.h>
#include <mkl_spblas.h>
#include <mkl_service.h>
#include <cmath>
#include <cstdio>
#include "matrixCSR.h"
#include "maxvol.h"
#include "nonlinearities.h"
#include "fullscheme.h"
#include "newton.h"
#include "reduction.h"
#include "constants.h"
#include "columbia.h"

int main(int argc, char ** argv)
{
	char basis[60];
	char nlbasis[60];
	char interpbasis[60];

	if (argc < 2)
	{
		printf("Input the basis storage file name tag.");
		return -2;
	};
	
	sprintf(basis, "bases/%sbasis.dat", argv[1]);
	sprintf(nlbasis, "bases/%sassist.dat", argv[1]);
	sprintf(interpbasis, "bases/%sinterp.dat", argv[1]);	

	builddata(GRIDsize, alphagrid, betagrid, shiftgrid, RANKmaxreduction, RANKnonlinearity, Snapshots, NLSnapshots, basis, nlbasis, interpbasis);
};

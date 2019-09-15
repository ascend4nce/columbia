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
	double * point, * refpoint;
	//double * nodes;
	int gridN;
	double shift, alpha, beta;
	
	double time;

	char basis[60];
	char bassist[60];
	char picturenamev[60];
	char picturenamex[60];
	char picturenamey[60];

	char bpicture1[60];
	char bpicture9[60];

	if (argc < 2)
	{
		printf("Please input the basis file name tag.\n");
		return -2;
	};

	sprintf(basis, "bases/%sbasis.dat", argv[1]);
	sprintf(bassist, "bases/%sassist.dat", argv[1]);

	sprintf(bpicture1, "baseprints/%snodes.dat", argv[1]);
	sprintf(bpicture9, "baseprints/%sdepnodes.dat", argv[1]);

	solverparams PRS(GRIDsize, Alpha, Beta, Shiftx, Shifty, Scaling);

	time = dsecnd();
	partsolverparams PPRS(RANKmaxreduction, RANKnonlinearity, bassist, &PRS);
	time = dsecnd() - time;
	//printf("Reduction preparations time: %f\n", time);

	PPRS.uploadbasis(basis);

	for(int i = 0; i < RANKreduction; i++)
	{
		sprintf(picturenamev, "baseprints/%svalue%d.dat", argv[1], i + 1);
		sprintf(picturenamex, "baseprints/%sfluxx%d.dat", argv[1], i + 1);
		sprintf(picturenamey, "baseprints/%sfluxy%d.dat", argv[1], i + 1);
		
		printsimple(&PRS, PPRS.V + i * PRS.varnum, picturenamev, picturenamex, picturenamey);
	};
	
	printdependencynodes(&PPRS, bpicture1, bpicture9);
	
	delete[] point;
	delete[] refpoint;
};

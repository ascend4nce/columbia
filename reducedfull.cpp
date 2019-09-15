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
	double * point, * smallpoint, * refpoint;
	double * nodes;
	int gridN;
	double shift, alpha, beta;
	FILE * results;
	FILE * known;

	results = fopen("reduction/fullresults.dat", "w");
	known = fopen("reduction/interpprint.dat", "w");	
	double time;

	char assist[60];
	char basis[60];
	char interp[60];

	if (argc < 2)
	{
		printf("Please input the basis file name tag.\n");
		return -2;
	};

	sprintf(assist, "bases/%sassist.dat", argv[1]);
	sprintf(basis, "bases/%sbasis.dat", argv[1]);
	sprintf(interp, "bases/%sinterp.dat", argv[1]);

	smallpoint = new double[RANKreduction];
	point = new double[6 * GRIDsize * GRIDsize + GRIDsize];
	refpoint = new double[6 * GRIDsize * GRIDsize + GRIDsize];
	
	double avglf = 0.0;
	double avglv = 0.0;
	double avgcf = 0.0;
	double avgcv = 0.0;
	double avgftmit = 0.0;
	double avgrtmit = 0.0;
	double vtavtm = 0.0;
	double loadtm = 0.0;

	srand((unsigned int) dsecnd());
	int Tests = betagrid * shiftgrid;
	for(int i = 0; i < Tests; i++)
	{
		double timefull, timesmall;
		solverparams * PRS = gridtoparams(GRIDsize,0, i % betagrid, i / betagrid, alphagrid, betagrid, shiftgrid);
		//solverparams * PRS = gridtoparams(GRIDsize,0, rand() % betagrid, rand() % shiftgrid, alphagrid, betagrid, shiftgrid);
	
		time = dsecnd();
		partsolverparams * PPRS = new partsolverparams(RANKmaxreduction, RANKnonlinearity, assist, PRS);
		PPRS->uploadinterpbasis(interp);
		time = dsecnd() - time;

		loadtm = time;

		PPRS->uploadbasis(basis);
		//fflush();
		time = dsecnd();     
		reducedtest(PPRS, RANKreduction, smallpoint);
		time = dsecnd() - time;
		timesmall = time;
		avgrtmit += timesmall / Tests;

		time = dsecnd();
		restoresolution(PPRS, smallpoint, point);
		time = dsecnd() - time;
		
		time = dsecnd();
		multitest(PRS, refpoint, GridNumber);
		time = dsecnd() - time;
		timefull = time;
		avgftmit += timefull / Tests;

		int fluxsize = 4 * GRIDsize * GRIDsize + GRIDsize;
		int vsize = 2 * GRIDsize * GRIDsize;
		int ione     = 1;
		double dmone = -1.0;

		double lferr, lverr;
		
		double cferr, cverr;

		if (i == Tests - 1)
		{
			printrs(PRS, "reduction/rs.dat");
			printreference(PRS, point, refpoint, "reduction/rd");
		};

		daxpy(&(PRS->varnum), &dmone, refpoint, &ione, point, &ione);
		//printf("Reduction relative flux error in 2 norm: %e.\n", dnrm2(&fluxsize, point, &ione) / dnrm2(&fluxsize, refpoint, &ione));
		lferr = dnrm2(&fluxsize, point, &ione) / dnrm2(&fluxsize, refpoint, &ione);
		//printf("Reduction relative value error in 2 norm: %e.\n", dnrm2(&vsize, point + fluxsize, &ione) / dnrm2(&vsize, refpoint + fluxsize, &ione));
		lverr = dnrm2(&vsize, point + fluxsize, &ione) / dnrm2(&vsize, refpoint + fluxsize, &ione);

		cferr = std::abs(point[idamax(&fluxsize, point, &ione) - 1]) / 
std::abs(refpoint[idamax(&fluxsize, refpoint, &ione) - 1]);

		cverr = std::abs(point[idamax(&vsize, point + fluxsize, &ione) - 1 + fluxsize]) / 
std::abs(refpoint[idamax(&vsize, refpoint + fluxsize, &ione) - 1 + fluxsize]);
		
		avglf += lferr / Tests;
		avglv += lverr / Tests;
		avgcf += cferr / Tests;
		avgcv += cverr / Tests;
		
		if (i == Tests - 1)
		{		
			time = dsecnd();
			PPRS->applyreductionshift(0);
			time = dsecnd() - time;
			vtavtm = time;
		};
		//PPRS->findoptimal(RANKreduction, refpoint, point);

		//daxpy(&(PRS->varnum), &dmone, refpoint, &ione, point, &ione);
		//printf("Optimal relative flux error in 2 norm: %e.\n", dnrm2(&fluxsize, point, &ione) / dnrm2(&fluxsize, refpoint, &ione));

		//printf("Optimal relative value error in 2 norm: %e.\n", dnrm2(&vsize, point + fluxsize, &ione) / dnrm2(&vsize, refpoint + fluxsize, &ione));
		if(i == Tests - 1)
		{
			for(int i = 0; i < PPRS->interpdata; i++)
			{
				fprintf(known,"%e %e %e\n", PPRS->betaexample[i], PPRS->shiftexample[i], 1.0);
			};
		};
		fprintf(results, "%e %e %e %e %e %e %f %f\n", PRS->beta, PRS->shifty, lferr, lverr, cferr, cverr, timefull, timesmall);

		delete PRS;
		delete PPRS;
	};

	if (printlevel >= 2)
	{
		printf("Reduction average results:\n");
		printf("l-2-error, flux: %e; l-2-error, values: %e; c-error, flux: %e; c-error, values: %e. \n", avglf, avglv, avgcf, avgcv);
		printf("Time: full %f; reduced %f; load %f; transform %f\n.", avgftmit, avgrtmit, loadtm, vtavtm);
	};

	fclose(results);
	fclose(known);

	delete[] point;
	delete[] refpoint;
	delete[] smallpoint;
};

#include <mkl.h>
#include <mkl_blas.h>
#include <cmath>
#include <cstdio>
#include "matrixCSR.h"
#include "maxvol.h"
#include "fullscheme.h"
#include "reduction.h"
#include "newton.h"
#include "constants.h"
#include "columbia.h"

void partsolverparams::findoptimal(int rk, double * vec, double * optapprox)
{
	char cN = 'N';
	char cT = 'T';

        double done  = 1.0;
        double dzero = 0.0;
        int ione     = 1;
        double * tmp;
        tmp = new double[rk];

	prs->vecshiftunroll(vec);
        dgemv(&cT, &(prs->varnum), &rk, &done, V, &(prs->varnum), vec, &ione, &dzero, tmp, &ione);
	prs->vecshiftrotate(vec);        
	dgemv(&cN, &(prs->varnum), &rk, &done, V, &(prs->varnum), tmp, &ione, &dzero, optapprox, &ione);
	prs->vecshiftrotate(optapprox);

        delete[] tmp;
};

void partsolverparams::testnl(int rk, int ntests)
{
	double * pot = new double[prs->varnum];
	double * smpt = new double[rk];
	double * nlmtx = new double[ntests * prs->varnum];
	double * svalues = new double[ntests];

	char cN = 'N';
	int ione = 1;
	double done = 1.0;
	double dzero = 0.0;
	int lwork = 64 * prs->varnum;
	double * work = new double[lwork];
	int info;

	for(int i = 0; i < ntests; i++)
	{
		for(int j = 0; j < rk; j++)
		{
			smpt[j] = 1.0e+5 * (-1.0 + 2.0 * (double) rand() / RAND_MAX);
		};
		
		dgemv(&cN, &(prs->varnum), &rk, &done, V, &(prs->varnum), smpt, &ione, &dzero, pot, &ione);
		prs->vecshiftrotate(pot);
		for(int j = 0; j < prs->varnum; j++)
		{
			nlmtx[i * prs->varnum + j] = oneresidnl(pot, (void *) prs, prs->shiftindex(j));
		};
		
	};

	dgesvd(&cN, &cN, &(prs->varnum), &ntests, nlmtx, &(prs->varnum), svalues, NULL, &(prs->varnum), NULL, &ntests, work, &lwork, &info);
	
	for(int i = 0; i < 1000; i++)
	{
		printf("%d: %e\n", i + 1, svalues[i] / svalues[0]);
	};
	delete[] svalues;
	delete[] work;
	
	delete[] nlmtx;
	delete[] pot;
	
	delete[] smpt;
	
};

void partsolverparams::uploadbasis(char * basis)
{
	FILE * fdb;
	int vnread, maxrank;
	double dummy;

	if (isLoaded == false)
	{
		V = new double[prs->varnum * rank];
		U = new double[prs->varnum * rank2];
	};

	fdb = fopen(basis, "r");	
	fscanf(fdb, "%d %d \n", &vnread, &maxrank);

	if ((vnread != prs->varnum) || (rank != maxrank))
	{
		printf("%d <> %d; %d <> %d \n\n", vnread, prs->varnum, rank, maxrank);
		printf("Invalid basis dimensionality\n");
		fclose(fdb);
		return;
	};

	for(int i = 0; i < prs->varnum; i++)
	{
		for(int j = 0; j < rank; j++)
		{
			fscanf(fdb, "%lf ", V + i + j * prs->varnum);
		};
		fscanf(fdb, "\n");
	};
	isLoaded = true;
	fclose(fdb);
};

void partsolverparams::uploadinterpbasis(char * basis)
{
	FILE * fdint;

	fdint = fopen(basis, "r");	
	fscanf(fdint, "%d\n", &interpdata);

	if (isInterpLoaded == false)
	{
		alphaexample = new double[interpdata];
		betaexample = new double[interpdata];
		shiftexample = new double[interpdata];

		interpmatrix = new double[interpdata * rank];
	};

	for(int i = 0; i < interpdata; i++)
	{
		fscanf(fdint, "%lf %lf %lf ", alphaexample + i, betaexample + i, shiftexample + i);
		//printf("%e %e %e \n", alphaexample[i], betaexample[i], shiftexample[i]);
		for(int j = 0; j < rank; j++)
		{
			fscanf(fdint, "%lf ", interpmatrix + i + j * interpdata);
		};	
		fscanf(fdint, "\n");
	};

	isInterpLoaded = true;
	fclose(fdint);
};

partsolverparams::partsolverparams(int rankset, int rank2set, char * assist, solverparams * prsset)
{
	int vn;
	int vnread, maxrank, baseread;

	FILE * fdassist;

	rank = rankset;
	rank2 = rank2set;	

	prs = prsset;
	vn = prs->varnum;
	
	isLoaded = false;
	isInterpLoaded = false;

        double dummy;
	int finalcheck;

	fdassist = fopen(assist, "r");

	//loading dimensions
	fscanf(fdassist, "%d %d %d %d\n\n", &vnread, &baseread, &maxrank, &rank9);

	if ((vnread != vn) || (rank2 != baseread) || (rank != maxrank))
	{
		fclose(fdassist);
		printf("Invalid basis parameters.\n");
		return;
	};

	current_rank = rank;

	ia1 = new int[rank2];
	ia9 = new int[rank9];
        
        VTAV  = new double[rank * rank];
        VTUSU = new double[rank * rank2];
        I9V   = new double[rank * rank9];

	Aepsex = new double[prs->n * rank];
	Aepsv = new double[prs->n * rank];

        //skipping svalues
        for(int i = 0; i < rank; i++)
        {
                fscanf(fdassist, "%lf ", &dummy);
        };
        fscanf(fdassist, "\n");
        for(int i = 0; i < rank2; i++)
        {
                fscanf(fdassist, "%lf ", &dummy);
        };
        fscanf(fdassist, "\n\n");

        //loading ia1, ia9
        for(int i = 0; i < rank2; i++)
        {
                fscanf(fdassist, "%d ", ia1 + i);
        };
        fscanf(fdassist, "\n");

        for(int i = 0; i < rank9; i++)
        {
                fscanf(fdassist, "%d ", ia9 + i);
        };
        fscanf(fdassist, "\n\n");

        //loading VTAV (part)
        for(int i = 0; i < rank; i++)
        {
                for(int j = 0; j < rank; j++)
                {
                        fscanf(fdassist,"%lf ", VTAV + i + j * rank);
                };
                fscanf(fdassist, "\n");
        };

	fscanf(fdassist, "\n");

	//loading shift support vectors for VTAV
	for(int i = 0; i < prs->n; i++)
	{
		for(int j = 0; j < rank; j++)
		{
			fscanf(fdassist, "%lf ", Aepsex + i + j * prs->n);
		};
		fscanf(fdassist, "\n");
	};
	fscanf(fdassist, "\n");

	for(int i = 0; i < prs->n; i++)
	{
		for(int j = 0; j < rank; j++)
		{
			fscanf(fdassist, "%lf ", Aepsv + i + j * prs->n);
		};
		fscanf(fdassist, "\n");
	};
	fscanf(fdassist, "\n");

	//loading VTUSU(part)
	for(int i = 0; i < rank; i++)
	{
		for(int j = 0; j < rank2; j++)
		{
			fscanf(fdassist, "%lf ", VTUSU + i + j * rank);
		};
		fscanf(fdassist, "\n");
	};
	fscanf(fdassist, "\n");

	//loading I9V(part)
	for(int i = 0; i < rank9; i++)
	{
		for(int j = 0; j < rank; j++)
		{
			fscanf(fdassist, "%lf ", I9V + i + j * rank9);
		};
		fscanf(fdassist, "\n");
	};
	fscanf(fdassist, "\n");

	//checking consistency
	fscanf(fdassist, "%d\n", &finalcheck);

	if (finalcheck != 0) printf("Basis file corrupted!\n");

	localshiftnum = 0;

	applyreductionshift(prs->shiftnum);

	fclose(fdassist);
};

void partsolverparams::applyreductionshift(int newshift)
{
	double positives;
	double negatives;

	positives = prs->scaling / prs->h;
	negatives = - prs->scaling / prs->h;
	
	for(int i = 0; i < prs->n; i++)
	{
		dger(&rank, &rank, &positives, Aepsex + i, &(prs->n), Aepsv + modL(i + localshiftnum, prs->n), &(prs->n), VTAV, &rank);
	};

	for(int i = 0; i < prs->n; i++)
	{
		dger(&rank, &rank, &negatives, Aepsex + i, &(prs->n), Aepsv + modL(i + newshift, prs->n), &(prs->n), VTAV, &rank);
	};

	for(int i = 0; i < prs->n; i++)
	{
		dger(&rank, &rank, &positives, Aepsv + modL(i + localshiftnum, prs->n), &(prs->n), Aepsex + i, &(prs->n), VTAV, &rank);
	};

	for(int i = 0; i < prs->n; i++)
	{
		dger(&rank, &rank, &negatives, Aepsv + modL(i + newshift, prs->n), &(prs->n), Aepsex + i, &(prs->n), VTAV, &rank);
	};

	localshiftnum = newshift;
};

void partsolverresidual(double * pt, double * ev, void * params)
{
	partsolverparams * PPRS;
	PPRS = (partsolverparams *) params;
	solverparams * PRS;
	PRS = PPRS->prs;

	char cN = 'N';
	char cT = 'T';
	double dzero = 0.0;
	double done = 1.0;
	double dmone = -1.0;
	int ione = 1;
	int izero = 0;

	int crk = PPRS->current_rank;
	int rk = PPRS->rank;
	int rk2 = PPRS->rank2;
	int rk9 = PPRS->rank9;
	int vn = PRS->varnum;
	int ctr;

	//linear residual part
	dgemv(&cN, &crk, &crk, &done, PPRS->VTAV, &rk, pt, &ione, &dzero, ev, &ione);

	//nonlinear residual part
	double * temp = new double[vn];
	double * partnl = new double[rk2];
	double * fullnl = new double[vn];

	double * test1 = new double[crk];
	double * test2 = new double[crk];

	for(int i = 0; i < rk9; i++)
	{
		temp[PRS->shiftindex(PPRS->ia9[i])] = ddot(&crk, PPRS->I9V + i, &rk9, pt, &ione); 
	};

	for(int i = 0; i < rk2; i++)
	{
		partnl[i] = oneresidnl(temp, (void *) PRS, PRS->shiftindex(PPRS->ia1[i]));
	};

	dgemv(&cN, &crk, &rk2, &done, PPRS->VTUSU, &rk, partnl, &ione, &done, ev, &ione);
		
	/*dgemv(&cN, &crk, &rk2, &done, PPRS->VTUSU, &rk, partnl, &ione, &dzero, test1, &ione);
	
	for(int i = 0; i < vn; i++)
	{
		temp[PRS->shiftindex(i)] = ddot(&crk, PPRS->V + i, &vn, pt, &ione); 
	};

	for(int i = 0; i < vn; i++)
	{
		fullnl[i] = oneresidnl(temp, (void *) PRS, PRS->shiftindex(i));
	};

	//dgemv(&cT, &vn, &crk, &done, PPRS->V, &vn, fullnl, &ione, &dzero, test2, &ione);

	//daxpy(&crk, &dmone, test2, &ione, test1, &ione);
	//printf("!%e!\n", dnrm2(&crk, test1, &ione) / dnrm2(&crk, test2, &ione));	
	*/
	//dgemv(&cN, &crk, &rk2, &done, PPRS->VTUSU, &rk, partnl, &ione, &dzero, ev, &ione);
	
	delete[] test1;
	delete[] test2;

	delete[] fullnl;
	delete[] temp;
	delete[] partnl;

};

void partsolverbadresidual(double * pt, double * ev, void * params)
{
	partsolverparams * PPRS;
	PPRS = (partsolverparams *) params;
	solverparams * PRS;
	PRS = PPRS->prs;

	char cN = 'N';
	char cT = 'T';
	double dzero = 0.0;
	double done = 1.0;
	double dmone = -1.0;
	int ione = 1;
	int izero = 0;

	int crk = PPRS->current_rank;
	int rk = PPRS->rank;
	int rk2 = PPRS->rank2;
	int rk9 = PPRS->rank9;
	int vn = PRS->varnum;
	int ctr;

	//linear residual part
	dgemv(&cN, &crk, &crk, &done, PPRS->VTAV, &rk, pt, &ione, &dzero, ev, &ione);

	//nonlinear residual part
	double * temp = new double[vn];
	double * fullnl = new double[vn];

	for(int i = 0; i < vn; i++)
	{
		temp[PRS->shiftindex(i)] = ddot(&crk, PPRS->V + i, &vn, pt, &ione); 
	};

	for(int i = 0; i < vn; i++)
	{
		fullnl[i] = oneresidnl(temp, (void *) PRS, PRS->shiftindex(i));
	};

	dgemv(&cT, &vn, &crk, &done, PPRS->V, &vn, fullnl, &ione, &done, ev, &ione);
	
	//dgemv(&cN, &crk, &rk2, &done, PPRS->VTUSU, &rk, partnl, &ione, &dzero, ev, &ione);
	
	delete[] temp;
	delete[] fullnl;

};



matrixCSR * partsolverderivnum(double * pt, void * params)
{
	partsolverparams * PPRS;
	PPRS = (partsolverparams *) params;
	solverparams * PRS = PPRS->prs;

	double * tempup, * tempdown, * temppoint;
	int ione = 1;
	int izero = 0;
	double done = 1.0;
	double dzero = 0.0;
	double dminusone = -1.0;

	int vn = PPRS->current_rank;
	double ash = 1.0e-4;

	matrixCSR * result;
	result    = new matrixCSR();
	result->n = vn;
	result->m = vn;
	result->a = new double[vn * vn];
	result->ia = new int[vn + 1];
	result->ja = new int[vn * vn];

	tempup = new double[vn];
	tempdown  = new double[vn];
	temppoint = new double[vn];

	dcopy(&vn, pt, &ione, temppoint, &ione);

	result->ia[vn] = vn * vn;
	for(int i = 0; i < vn; i++)
	{
		result->ia[i] = i * vn;
		for(int j = 0; j < vn; j++)
		{
			result->ja[i * vn + j] = j;
		};
	};

	for(int j = 0; j < vn; j++)
	{
		double add;
		double addinv;

		add = ash * dnrm2(&vn, temppoint, &ione);
		add = ash;
		addinv = 1.0 / add;

		temppoint[j] += add  / 2.0;
		partsolverresidual(temppoint, tempup, params);
		//partsolverbadresidual(temppoint, tempup, params);
		temppoint[j] -= add;

		partsolverresidual(temppoint, tempdown, params);
		//partsolverbadresidual(temppoint, tempdown, params);
		temppoint[j] += add / 2.0;

		daxpy(&vn, &dminusone, tempdown, &ione, tempup, &ione);
		dscal(&vn, &addinv, tempup, &ione);

		dcopy(&vn, tempup, &ione, result->a + j, &vn);
		//daxpy(&vn, pt + j, tempup, &ione, a + j , &vn);
	};

	for(int i = 0; i < vn + 1; i++)
	{
		result->ia[i]+=1;
	};
	for(int i = 0;i < vn * vn; i++)
	{
		result->ja[i]+=1;
	};

	delete[] temppoint;
	delete[] tempup;
	delete[] tempdown;

	return result;
};

matrixCSR * partsolverderivCSR(double * pt, void * params)
{
	partsolverparams * PPRS;
	PPRS = (partsolverparams *) params;
	solverparams * PRS;
	PRS = PPRS->prs;

	char cN = 'N';
	char cT = 'T';
	double dzero = 0.0;
	double done = 1.0;
	double dmone = -1.0;
	int ione = 1;
	int izero = 0;

	int rk = PPRS->rank;
	int crk = PPRS->current_rank;
	int rk2 = PPRS->rank2;
	int rk9 = PPRS->rank9;
	int vn = PRS->varnum;

	double * matrix = new double[crk * crk];
	double * minijacob = new double[rk2 * crk];
	double * temp = new double[vn];
	int * shiftedia1 = new int[rk2];

	for(int i = 0; i < crk; i++)
	{
		dcopy(&crk, PPRS->VTAV + i, &rk, matrix + i, &crk);
	};

	for(int i = 0; i < rk9; i++)
	{
		temp[PRS->shiftindex(PPRS->ia9[i])] = ddot(&crk, PPRS->I9V + i, &rk9, pt, &ione); 
	};
	
	matrixCSR * smallderiv;
	for(int i = 0; i < rk2; i++)
	{
		shiftedia1[i] = PRS->shiftindex(PPRS->ia1[i]);
	};
	smallderiv = partderivCSR(temp, PRS, rk2, shiftedia1);

	for(int j = 0; j < crk; j++)
	{
		for(int i = 0; i < rk9; i++)
		{
			temp[PRS->shiftindex(PPRS->ia9[i])] = PPRS->I9V[i + j * rk9];
		};
		smallderiv->matvec(temp, minijacob + rk2 * j);
	};

	delete smallderiv;

	dgemm(&cN, &cN, &crk, &crk, &rk2, &done, PPRS->VTUSU, &rk, minijacob, &rk2, &done, matrix, &crk);
	//dgemm(&cN, &cN, &crk, &crk, &rk2, &done, PPRS->VTUSU, &rk, minijacob, &rk2, &dzero, matrix, &crk);


	matrixCSR * result;
	result    = new matrixCSR();
	result->n = crk;
	result->m = crk;
	result->a = new double[crk * crk];
        result->ia = new int[crk + 1];
        result->ja = new int[crk * crk];

        result->ia[crk] = crk * crk;
        for(int i = 0; i < crk; i++)
        {
                result->ia[i] = i * crk;
                for(int j = 0; j < crk; j++)
                {
                        result->ja[i * crk + j] = j;
			result->a[i * crk + j] = matrix[i + crk * j];
                };
        };

	for(int i = 0; i < crk + 1; i++)
	{
		result->ia[i]+=1;
	};
	for(int i = 0;i < crk * crk; i++)
	{
		result->ja[i]+=1;
	};

	delete[] shiftedia1;
	delete[] temp;
	delete[] minijacob;
	delete[] matrix;
	
	return result;
};


partsolverparams::~partsolverparams()
{
	if (isLoaded)
	{
		delete[] V;
	};
	if (isInterpLoaded)
	{
		delete[] alphaexample;
		delete[] betaexample;
		delete[] shiftexample;
		delete[] interpmatrix;
	};
	delete[] ia9;
	delete[] ia1;
	delete[] VTAV;
	delete[] VTUSU;
	delete[] I9V;
	delete[] Aepsex;
	delete[] Aepsv;
};

void reducedtest(partsolverparams * PPRS, int rk, double * smallpoint)
{	
	solverparams * PRS = PPRS->prs;

	PPRS->current_rank = rk;

	double * smpoint2;	
	double zeroresidual;
	double obtresidual;

	double tlr;
	int ione = 1;

	smpoint2 = new double[PPRS->current_rank];

	for(int i = 0; i < PPRS->current_rank; i++)
	{
		smallpoint[i] = 1.0e-4 / PRS->scaling;
	};

	partsolverresidual(smallpoint, smpoint2, PPRS);
	zeroresidual = dnrm2(&(PPRS->current_rank), smpoint2, &ione);

	if (PPRS->isInterpLoaded)
	{
		int optnum = -1;
		double distmin;

		for(int j = 0; j < PPRS->interpdata; j++)
		{
			double dist = 0;	
			dist += std::abs(PRS->alpha - PPRS->alphaexample[j]);
			dist += std::abs(PRS->beta - PPRS->betaexample[j]) * PRS->scaling;
			dist+= std::abs(PRS->shifty - PPRS->shiftexample[j]);

			if ((optnum == -1) || (dist < distmin))
			{
				optnum = j;
				distmin = dist;
			};
		};

		for(int i = 0; i < PPRS->current_rank; i++)
		{
			smallpoint[i] = PPRS->interpmatrix[optnum + i * PPRS->interpdata];

		};
	};

	partsolverresidual(smallpoint, smpoint2, PPRS);
	obtresidual = dnrm2(&(PPRS->current_rank), smpoint2, &ione);

	if (reducednewtonmethod == 2)
	{
		if (reducedderivatives == 1)
		{
			newtonmark2(PPRS->current_rank, smallpoint, partsolverresidual, partsolverderivCSR, PPRS, Tolerance * zeroresidual / obtresidual, MAXnewtoniterations, 10.0 / PRS->scaling,  LStolerance, METHOD_SVD);
		}
		else
		{
			newtonmark2(PPRS->current_rank, smallpoint, partsolverresidual, partsolverderivnum, PPRS, Tolerance * zeroresidual / obtresidual, MAXnewtoniterations, 10.0 / PRS->scaling,  LStolerance, METHOD_SVD);
		};
	}
	else
	{
		if (reducedderivatives == 1)
		{
			newton(PPRS->current_rank, smallpoint, partsolverresidual, partsolverderivCSR, PPRS, Tolerance * zeroresidual / obtresidual, MAXnewtoniterations, 10.0 / PRS->scaling,  LStolerance, METHOD_SVD, METHOD_COVSVD);
		}
		else
		{
			newton(PPRS->current_rank, smallpoint, partsolverresidual, partsolverderivnum, PPRS, Tolerance * zeroresidual / obtresidual, MAXnewtoniterations, 10.0 / PRS->scaling,  LStolerance, METHOD_SVD, METHOD_COVSVD);
		};
	};

	//printf("%e %e really?\n\n", zeroresidual, obtresidual);
	//newtonmark2(PPRS->current_rank, smallpoint, partsolverresidual, partsolverderivCSR, PPRS, Tolerance * zeroresidual / obtresidual, MAXnewtoniterations, 10.0 / PRS->scaling, LStolerance, METHOD_SVD);
	//newton(PPRS->current_rank, smallpoint, partsolverresidual, partsolverderivCSR, PPRS, Tolerance * zeroresidual / obtresidual, MAXnewtoniterations, 10.0 / PRS->scaling, LStolerance, METHOD_SVD, METHOD_COVSVD);
	//newton(PPRS->current_rank, smallpoint, partsolverbadresidual, partsolverderivnum, PPRS, Tolerance * zeroresidual / obtresidual, MAXnewtoniterations, 10.0 / PRS->scaling, LStolerance, METHOD_SVD, METHOD_COVSVD);


	for(int i = 0; i < PPRS->current_rank; i++)
	{
		if (printlevel >= 2) printf("(%e) ", smallpoint[i]);
	};
	if(printlevel >= 2) printf("\n");

	delete[] smpoint2;
};

//was a failed experiment to iteratively increase ranks. Experiment gave no results, do not use that function.
void reducedtestmark2(partsolverparams * PPRS, int rk,double * smallpoint)
{
	solverparams * PRS = PPRS->prs;
	smallpoint[0] = 1.0e+6;

	for(int i = 1; i < rk; i++)
	{
		PPRS->current_rank = i + 1;
		smallpoint[i] = 1.0e-4;

		newtonmark2(PPRS->current_rank, smallpoint, partsolverresidual, partsolverderivCSR, PPRS, Tolerance, MAXnewtoniterations, 10.0 / PRS->scaling, LStolerance, METHOD_SVD);
		//newton(PPRS->current_rank, smallpoint, partsolverresidual, partsolverderivCSR, PPRS, Tolerance, MAXnewtoniterations, 10.0 / PRS->scaling, LStolerance, METHOD_SVD, METHOD_COVSVD);

		for(int i = 0; i < PPRS->current_rank; i++)
		{
			printf("(%e) ", smallpoint[i]);
		};
		printf("\n");

	};
	
};

void restoresolution(partsolverparams * PPRS, double * smallpoint, double * point)
{
	solverparams * PRS = PPRS->prs;

	int izero    = 0, ione = 1;
	double dzero = 0.0;
	double done = 1.0;
	char cN = 'N';

	if (!(PPRS->isLoaded)) return;	

	dgemv(&cN, &(PRS->varnum), &(PPRS->current_rank), &done, PPRS->V, &(PRS->varnum), smallpoint, &ione, &dzero, point, &ione);
	PRS->vecshiftrotate(point);

};


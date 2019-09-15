#include <mkl.h>
#include <mkl_rci.h>
#include <mkl_blas.h>
#include <mkl_spblas.h>
#include <mkl_service.h>
#include <cmath>
#include <cstdio>
#include "matrixCSR.h"
#include "fullscheme.h"
#include "reduction.h"
#include "newton.h"
#include "maxvol.h"
#include "constants.h"
#include "columbia.h"

void simpletest(solverparams * PRS, double * point)
{	
	int izero    = 0, ione = 1;
	double dzero = 1e-12;
	
	dcopy(&(PRS->varnum), &dzero, &izero, point, &ione);
	
	if (fullnewtonmethod == 2)
	{
		if (fullderivatives == 1)
		{
			newtonmark2(PRS->varnum, point, fullsolverresidual, fullsolverderivCSR, PRS, Tolerance, MAXnewtoniterations, 10.0 / PRS->scaling,  LStolerance, METHOD_PRECGMRES);
		}
		else
		{
			newtonmark2(PRS->varnum, point, fullsolverresidual, fullsolverderivnum, PRS, Tolerance, MAXnewtoniterations, 10.0 / PRS->scaling,  LStolerance, METHOD_PRECGMRES);
		};
	}
	else
	{
		if (fullderivatives == 1)
		{
			newton(PRS->varnum, point, fullsolverresidual, fullsolverderivCSR, PRS, Tolerance, MAXnewtoniterations, 10.0 / PRS->scaling,  LStolerance, METHOD_PRECGMRES, METHOD_PRECCOVCG);
		}
		else
		{
			newton(PRS->varnum, point, fullsolverresidual, fullsolverderivnum, PRS, Tolerance, MAXnewtoniterations, 10.0 / PRS->scaling,  LStolerance, METHOD_PRECGMRES, METHOD_PRECCOVCG);
		};
	};
	
};

void pseudoinverse(int sz, double * mtx)
{
	int lwork = 16 * sz;
	double * work = new double[lwork];
	double * sv = new double[sz];
	double * u = new double[sz * sz];
	double * vt = new double[sz * sz];

	char cA = 'A';
	char cN = 'N';
	char cT = 'T';
	int ione = 1;
	int izero = 0;
	double done = 1.0;
	double dzero = 0.0;

	double scal;

	int info;
	int sz2 = sz * sz;
	
	dgesvd(&cA, &cA, &sz, &sz, mtx, &sz, sv, u, &sz, vt, &sz, work, &lwork, &info);
//	dcopy(&sz2, &dzero, &izero, mtx, &ione);

	for(int i = 0; i < sz; i++)
	{
		scal = 1.0 / sv[i];
		//printf("%d: %e\n", i, sv[i]);
//		dger(&sz, &sz, &scal, vt + i, &sz, u + i * sz, &ione, mtx, &sz);
                dscal(&sz, &scal, u + i * sz, &ione);
	};
        dgemm(&cT, &cT, &sz, &sz, &sz, &done, vt, &sz, u, &sz, &dzero, mtx, &sz);

	delete[] u;
	delete[] vt;
	delete[] sv;
	delete[] work;
};

void multitest(solverparams * PRSinit, double * point,  int doubling)
{
	int izero = 0, ione = 1;
	double dzero = 1e-12;

	double * pointsmall;
	double * pointlarge;

	double stepresidual;

	solverparams * PRSsmall;
	solverparams * PRSlarge;

	int finishn, startn = PRSinit->n;
	for(int i = 0; i < doubling - 1; i++)
	{
		startn /= 2.0;
	};

	finishn = startn;
	for(int i = 0; i < doubling - 1; i++)
	{
		finishn *= 2.0;
	};

	if (finishn != PRSinit->n) 
	{
		printf("Columbia: wrong grid hierarchy parameters\n");
		return;
	};
	
	PRSsmall = new solverparams(startn, PRSinit->alpha, PRSinit->beta, PRSinit->shiftx, PRSinit->shifty, PRSinit->scaling);
	pointsmall = new double[PRSsmall->varnum];

	dcopy(&(PRSsmall->varnum), &dzero, &izero, pointsmall, &ione);
	
	for(int i = 0; i < doubling; i++)
	{
		//calculating required residual
		double * test = new double[PRSsmall->varnum];
		double * testev = new double[PRSsmall->varnum];

		double normin, normzero;
		double gridtime;

		dcopy(&(PRSsmall->varnum), &dzero, &izero, test, &ione);
		fullsolverresidual(test, testev, PRSsmall);
		normzero = dnrm2(&(PRSsmall->varnum), testev, &ione);

		fullsolverresidual(pointsmall, testev, PRSsmall);
		normin = dnrm2(&(PRSsmall->varnum), testev, &ione);

		if (printlevel >= 2) printf("%e %e \n", normzero, normin);
		if (printlevel >= 1) printf("\n\n*** GRID %d x %d, starting residual %e \n\n", PRSsmall->n, PRSsmall->n * 2, normin/normzero);

		gridtime = dsecnd();

		if (fullnewtonmethod == 2)
		{
			if (fullderivatives == 1)
			{
				newtonmark2(PRSsmall->varnum, pointsmall, fullsolverresidual, fullsolverderivCSR, PRSsmall, Tolerance * normzero / normin, MAXnewtoniterations, 10.0 / PRSinit->scaling,  LStolerance, METHOD_PRECGMRES);
			}
			else
			{
				newtonmark2(PRSsmall->varnum, pointsmall, fullsolverresidual, fullsolverderivnum, PRSsmall, Tolerance * normzero / normin, MAXnewtoniterations, 10.0 / PRSinit->scaling,  LStolerance, METHOD_PRECGMRES);
			};
		}
		else
		{
			if (fullderivatives == 1)
			{
				newton(PRSsmall->varnum, pointsmall, fullsolverresidual, fullsolverderivCSR, PRSsmall, Tolerance * normzero / normin, MAXnewtoniterations, 10.0 / PRSinit->scaling,  LStolerance, METHOD_PRECGMRES, METHOD_PRECCOVCG);
			}
			else
			{
				newton(PRSsmall->varnum, pointsmall, fullsolverresidual, fullsolverderivnum, PRSsmall, Tolerance * normzero / normin, MAXnewtoniterations, 10.0 / PRSinit->scaling,  LStolerance, METHOD_PRECGMRES, METHOD_PRECCOVCG);
			};
		};

		//newtonmark2(PRSsmall->varnum, pointsmall, fullsolverresidual, fullsolverderivCSR, PRSsmall, Tolerance * normzero / normin, MAXnewtoniterations, 10.0 / PRSinit->scaling, LStolerance, METHOD_PRECGMRES);
		//newtonmark2(PRSsmall->varnum, pointsmall, fullsolverresidual, fullsolverderivCSR, PRSsmall, Tolerance * normzero / normin, MAXnewtoniterations, 10.0 / PRSinit->scaling, LStolerance, METHOD_SVD);


		if (i < doubling - 1)
		{		
			PRSlarge = new solverparams(PRSsmall->n * 2, PRSinit->alpha, PRSinit->beta, PRSinit->shiftx, PRSinit->shifty, PRSinit->scaling);
			pointlarge = new double[PRSlarge->varnum];
			int gn = PRSsmall->n;
			for(int i = 0; i < 2 * gn; i++)
			{
				for(int j = 0; j < gn; j++)
				{
					pointlarge[PRSlarge->pnum(2 * i, 2 * j)] = pointsmall[PRSsmall->pnum(i,j)];
					pointlarge[PRSlarge->pnum(2 * i + 1, 2 * j)] = pointsmall[PRSsmall->pnum(i,j)];
					pointlarge[PRSlarge->pnum(2 * i, 2 * j + 1)] = pointsmall[PRSsmall->pnum(i,j)];
					pointlarge[PRSlarge->pnum(2 * i + 1, 2 * j + 1)] = pointsmall[PRSsmall->pnum(i,j)];
				};
			};

			for(int j = 0; j < gn; j++)
			{
				for(int i = 0; i < 2 * gn; i++)
				{
					double left, right;
					left = pointsmall[PRSsmall->exnum(i,j)];
					right = pointsmall[PRSsmall->exnum(i + 1, j)];

					pointlarge[PRSlarge->exnum(2 * i, 2 * j)] = left;
					pointlarge[PRSlarge->exnum(2 * i + 1, 2 * j)] = (left + right) / 2.0;
					pointlarge[PRSlarge->exnum(2 * i, 2 * j + 1)] = left;
					pointlarge[PRSlarge->exnum(2 * i + 1, 2 * j + 1)] = (left + right) / 2.0;
				};

				pointlarge[PRSlarge->exnum(4 * gn, 2 * j)] = pointsmall[PRSsmall->eynum(2 * gn, j)];
				pointlarge[PRSlarge->exnum(4 * gn, 2 * j + 1)] = pointsmall[PRSsmall->eynum(2 * gn, j)];
			};

			for(int i = 0; i < 2 * gn; i++)
			{
				for(int j = 0; j < gn; j++)
				{
					double top, btm;
					btm = pointsmall[PRSsmall->eynum(i,j)];
					top = pointsmall[PRSsmall->eynum(i, j + 1)];

					pointlarge[PRSlarge->eynum(2 * i, 2 * j)] = btm;
					pointlarge[PRSlarge->eynum(2 * i + 1, 2 * j)] = btm;
					pointlarge[PRSlarge->eynum(2 * i, 2 * j + 1)] = (top + btm) / 2.0;
					pointlarge[PRSlarge->eynum(2 * i + 1, 2 * j + 1)] = (top + btm) / 2.0;
				};
				pointlarge[PRSlarge->eynum(2 * i, 2 * gn)] = pointsmall[PRSsmall->eynum(i, gn)];
				pointlarge[PRSlarge->eynum(2 * i + 1, 2 * gn)] = pointsmall[PRSsmall->eynum(i, gn)];
			};

			delete[] pointsmall;
			pointsmall = pointlarge;

			delete PRSsmall;
			PRSsmall = PRSlarge;

		};

		gridtime = dsecnd() - gridtime;

		int seconds, minutes, hours;

		hours = gridtime / 3600;
		minutes = gridtime / 60 - hours * 60;
		seconds = gridtime - minutes * 60 - hours * 3600;
	
		if (printlevel >= 1)
		{
			printf("\n\n Grid time:");
			if (hours > 0)
			{
				printf(" %d hours,", hours);
			};
			if (minutes > 0)
			{
				printf(" %d minutes,", minutes);
			};
			printf(" %d seconds.\n\n", seconds);
		};

		delete[] test;
		delete[] testev;
	};	

	delete PRSsmall;

	dcopy(&(PRSinit->varnum), pointsmall, &ione, point, &ione);
	
	delete[] pointsmall;
};

solverparams * gridtoparams(int gridN,int a, int b, int s, int as, int bs, int ss)
{
	double alpha, beta, shift;
	alpha = 1.6 * ((double) a / as);
	beta = 3.0e+4 * pow(100.0, (double) b / bs);
	shift = 1.0 * ((double) s / ss);
	solverparams * PRS = new solverparams(gridN, alpha, beta, Shiftx, shift, Scaling);
	return PRS;
};

void setsparsegrid(int seed, int n, int * a, int * b, int * s, int as, int bs, int ss)
{
	srand(seed);

	for(int i = 0; i < n; i++)
	{
		bool ok;
		do
		{
			ok = true;

			a[i] = rand() % as;
			b[i] = rand() % bs;
			s[i] = rand() % ss;
		
			for(int j = 0; j < i; j++)
			{
				if (a[i] == a[j] && b[i] == b[j] && s[i] == s[j])
				{
					ok = false;
				};
			};
		}
		while (!ok);
	};
};


void builddata(int gridN, int as, int bs, int ss, int rankmax, int basenum, int snapshots, int nlsnapshots, char * basis, char * assist, char * interp)
{
	FILE * fdb, * fdassist, * fdint;

	double * vmatrix;
	double * nlmatrix;

	int lwork;
	double * work;
	int info;
	char cO = 'O';
	char cN = 'N';
        char cT = 'T';
	char cS = 'S';
	double * svalues;
	double * snlvalues;

	int izero    = 0;
	int ione     = 1;
	double dzero = 0.0;
	double done  = 1.0;

	int * ia1, * ia9;
        int * ipiv;
        int rank9;

        double * VTAV, * VTUSU;
	double * intsvd;

	double * smpt;

        double * temp, * temp2, * temp3;
	
	VSLStreamStatePtr vstream;

	int * a, * b, * s;

	a = new int[snapshots + 1];
	b = new int[snapshots + 1];
	s = new int[snapshots + 1];

	svalues = new double[snapshots];
	snlvalues = new double[nlsnapshots];
	
	solverparams PRS(gridN, 0.0, 0.0, Shiftx, 0.0, Scaling);

	lwork = PRS.varnum * 16;
	work = new double[lwork];

	vmatrix = new double[PRS.varnum * snapshots];
	nlmatrix = new double[PRS.varnum * nlsnapshots];

	intsvd = new double[snapshots * snapshots];

	setsparsegrid(776, snapshots + 1, a, b, s, as, bs, ss);

	if (printlevel >= 2) printf("Not used solution: %d %d %d.\n", a[snapshots], b[snapshots], s[snapshots]);

	for(int j = 0; j < snapshots; j++)
	{
		solverparams * PRSgrid;

		PRSgrid = gridtoparams(PRS.n,a[j],b[j],s[j],as,bs,ss);
	
		if (printlevel >= 1) printf("%d: alpha %f, beta %f, shifty %f \n",j , PRSgrid->alpha, PRSgrid->beta, PRSgrid->shifty);

		multitest(PRSgrid, vmatrix + j * PRS.varnum, GridNumber);

		PRSgrid->vecshiftunroll(vmatrix + j * PRS.varnum);
		
		delete PRSgrid;
	};

	dgesvd(&cO, &cS, &(PRS.varnum), &snapshots, vmatrix, &(PRS.varnum), svalues, NULL, &(PRS.varnum), intsvd, &snapshots, work, &lwork, &info);
	if (info != 0) printf("Builder: dgesvd failed!\n");

	smpt = new double[rankmax];
	temp = new double[PRS.varnum];

	vslNewStream(&vstream, VSL_BRNG_MT2203, 776);
	
	for(int i = 0; i < nlsnapshots; i++)
	{
		solverparams * PRSNL;

		PRSNL = gridtoparams(PRS.n,rand() % as,rand() % bs, rand() % ss,as,bs,ss);

		for(int j = 0; j < rankmax; j++)
		{
			vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, vstream, 1, smpt + j, 0, svalues[j]);

			//vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, vstream, 1, smpt + j, 0, 1.0e+5);
		};

		dgemv(&cN, &(PRS.varnum), &rankmax, &done, vmatrix, &(PRS.varnum), smpt, &ione, &dzero, temp, &ione);

		PRSNL->vecshiftrotate(temp);

		for(int j = 0; j < PRS.varnum; j++)
		{
			nlmatrix[i * PRS.varnum + j] = oneresidnl(temp, (void *) PRSNL, PRSNL->shiftindex(j));
		};

		delete PRSNL;
	};

	vslDeleteStream(&vstream);

	delete[] temp;
	delete[] smpt;

	dgesvd(&cO, &cN, &(PRS.varnum), &nlsnapshots, nlmatrix, &(PRS.varnum), snlvalues, NULL, &(PRS.varnum), NULL, &nlsnapshots, work, &lwork, &info);
	if (info != 0) printf("Builder: dgesvd failed!\n");

	fdb = fopen(basis, "w");

	fprintf(fdb, "%d %d\n", PRS.varnum, rankmax);	
	
	for(int i = 0; i < PRS.varnum; i++)
	{
		for(int j = 0; j < rankmax; j++)
		{
			fprintf(fdb,"%e ", vmatrix[i + j * PRS.varnum]);
		};
		fprintf(fdb, "\n");
	};

        fclose(fdb);

	fdint = fopen(interp, "w");

	fprintf(fdint, "%d\n", snapshots);

	for(int j = 0; j < snapshots; j++)
	{
		solverparams * PRSgrid;

		PRSgrid = gridtoparams(PRS.n,a[j],b[j],s[j],as,bs,ss);

		fprintf(fdint,"%e %e %e ", PRSgrid->alpha, PRSgrid->beta, PRSgrid->shifty);

		for(int k = 0; k < rankmax; k++)
		{
			fprintf(fdint,"%e ", intsvd[k + j * snapshots] * svalues[k]);
		};

		fprintf(fdint, "\n");
		delete PRSgrid;
	};
	fclose(fdint);
	
        ia1 = new int[PRS.varnum];
        ia9 = new int[PRS.varnum];

	for(int i = 0; i < PRS.varnum; i++)
	{
		ia1[i] = i;
	};

	for(int k = 0; k < basenum; k++)
	{
		int tmp;
		int pivotx, pivoty, pivot;

		pivotx = rand() % (2 * PRS.n + 1);
		pivoty = rand() % PRS.n;

		pivot  = PRS.exnum(pivotx, pivoty);

                int sch;

                for(sch = 0; sch < PRS.nlnum; sch++)
                {
                        if (ia1[sch] == pivot)
                        {
                                break;
                        };
                };
                ia1[sch] = ia1[k];
                ia1[k]   = pivot;
	};

        maxvol2(basenum, PRS.varnum, nlmatrix, PRS.varnum, ia1);
	
        rank9 = 0;
        
        for(int i = 0; i < PRS.varnum; i++)
        {
                ia9[i] = 0;
        };

        for(int j = 0; j < PRS.n; j++)
        {
                ia9[PRS.exnum(PRS.halfnum, j)] = 1.0;
                ia9[PRS.exnum(PRS.halfnum + 1, j)] = 1.0;
                ia9[PRS.eynum(PRS.halfnum, j)] = 1.0;
                rank9 += 3;
        };

        matrixCSR * ootka;
        ootka = partderivCSR(vmatrix, (void *) &PRS, basenum, ia1);

	//ootka->print();
        for(int i = 0; i < ootka->ia[basenum] - 1; i++)
        {
                if (ia9[ootka->ja[i] - 1] == 0) rank9++;
                ia9[ootka->ja[i] - 1] = 1;
        };

        delete ootka;

        VTAV  = new double[rankmax * rankmax];
        VTUSU = new double[rankmax * basenum];

        temp  = new double[PRS.varnum];
        temp2 = new double[basenum * basenum];
        temp3 = new double[rankmax * basenum];

        for(int i = 0; i < basenum; i++)
        {
                dcopy(&basenum, nlmatrix + ia1[i], &(PRS.varnum), temp2 + i, &basenum);
        };

        ipiv = new int[basenum];
	/*
        //do svd here instead?>
        dgetrf(&basenum, &basenum, temp2, &basenum, ipiv, &info);
        if (info != 0) printf("dgetrf failed!\n");
        lwork = rankmax * basenum;
        dgetri(&basenum, temp2, &basenum, ipiv, temp3, &lwork, &info);
        if (info != 0) printf("dgetri failed!\n");
	*/
	pseudoinverse(basenum, temp2);

        dgemm(&cT, &cN, &rankmax, &basenum, &(PRS.varnum), &done, vmatrix, &(PRS.varnum), nlmatrix, &(PRS.varnum), &dzero, temp3, &rankmax);
        dgemm(&cN, &cN, &rankmax, &basenum, &basenum, &done, temp3, &rankmax, temp2, &basenum, &dzero, VTUSU, &rankmax);

        for(int j = 0; j < rankmax; j++)
        {
                for(int k = 0; k < PRS.varnum; k++)
                {
                        temp[k] = oneresidlin(vmatrix + j * PRS.varnum, (void *) &PRS, k);
                };
                dgemv(&cT, &(PRS.varnum), &rankmax, &done, vmatrix, &(PRS.varnum), temp, &ione, &dzero, VTAV + j * rankmax, &ione);
        };

	fdassist = fopen(assist, "w");

	//saving dimensions
	fprintf(fdassist, "%d %d %d %d\n\n", PRS.varnum, basenum, rankmax, rank9);

	//saving simgular values of V and U
	for(int i = 0; i < rankmax; i++)
	{
		fprintf(fdassist, "%e ", svalues[i]);
	};
	fprintf(fdassist, "\n");

	for(int i = 0; i < basenum; i++)
	{
		fprintf(fdassist, "%e ", snlvalues[i]);
	};
	fprintf(fdassist, "\n");

	fprintf(fdassist, "\n");

	//saving maxvol results
	for(int i = 0; i < basenum; i++)
	{
		fprintf(fdassist, "%d ", ia1[i]);
	};
	fprintf(fdassist, "\n");

	int count =  0;
	for(int i = 0; i < PRS.varnum; i++)
	{
		if (ia9[i] == 1)
		{
			count++;
			fprintf(fdassist, "%d ", i);
		};
	};
	if (count != rank9) 
	{
		printf("Builder: rank9 mismatch\n");
	};
	fprintf(fdassist, "\n\n");

	//svaing VTAV
	for(int i = 0; i < rankmax; i++)
	{
		for(int j = 0; j < rankmax; j++)
		{
			fprintf(fdassist, "%e ", VTAV[i + j * rankmax]);
		};
		fprintf(fdassist, "\n");
	};

	fprintf(fdassist, "\n");
	
	//saving shift support vectors for VTAV
	for(int i = 0; i < PRS.n; i++)
	{
		for(int j = 0; j < rankmax; j++)
		{
			fprintf(fdassist, "%e ", vmatrix[PRS.exnum(PRS.halfnum, i) + j * PRS.varnum]);
		};
		fprintf(fdassist, "\n");
	};
	fprintf(fdassist, "\n");

	for(int i = 0; i < PRS.n; i++)
	{
		for(int j = 0; j < rankmax; j++)
		{
			fprintf(fdassist, "%e ", vmatrix[PRS.pnum(PRS.halfnum, i) + j * PRS.varnum]);
		};
		fprintf(fdassist, "\n");
	};
	fprintf(fdassist, "\n");
	
	//saving VTUSU
	for(int i = 0; i < rankmax; i++)
	{
		for(int j = 0; j < basenum; j++)
		{
			fprintf(fdassist, "%e ", VTUSU[i + j * rankmax]);
		};
		fprintf(fdassist, "\n");
	};

	fprintf(fdassist, "\n");

	//saving I9V
	for(int i = 0; i < PRS.varnum; i++)
	{
		if (ia9[i] == 1)
		{
			for(int j = 0; j < rankmax; j++)
			{
				fprintf(fdassist, "%e ", vmatrix[i + j * PRS.varnum]);
			};
			fprintf(fdassist, "\n");
		};
	};

	fprintf(fdassist, "\n0\n");
	
	fclose(fdassist);

	delete[] a;
	delete[] b;
	delete[] s;

        delete[] ia1;
        delete[] ia9;
        delete[] ipiv;

        delete[] VTAV;
        delete[] VTUSU;
	delete[] intsvd;

        delete[] temp;
        delete[] temp2;
        delete[] temp3;
		
	delete[] work;
	delete[] svalues;
	delete[] snlvalues;
	delete[] vmatrix;
	delete[] nlmatrix;

};

void printdependencynodes(partsolverparams * PPRS, char * where1, char * where9)
{
	solverparams * PRS = PPRS->prs;
	FILE * fd1, * fd9;
	fd1 = fopen(where1, "w");
	fd9 = fopen(where9, "w");

	double x, y;

	double * pt1 = new double[PRS->varnum];
	double * pt9 = new double[PRS->varnum];

	int gridN = PRS->n;
	
	for(int i = 0; i < PRS->varnum; i++)
	{
		pt1[i] = 1.0;
		pt9[i] = 1.0;
	};

	for(int i = 0; i < PPRS->rank2; i++)
	{
		pt1[PPRS->ia1[i]] = 0.0;
	};

	for(int i = 0; i < PPRS->rank9; i++)
	{
		pt9[PPRS->ia9[i]] = 0.0;
	};

	//values
	for(int i = 0; i < 2 * gridN; i++)
	{
		for(int j = 0; j < gridN; j++)
		{
			x = PRS->h / 2.0 + i * PRS->h;
			y = PRS->h / 2.0 + j * PRS->h;	

			fprintf(fd1, "%f %f %e\n", x, y, pt1[PRS->pnum(i,j)]);
			fprintf(fd9, "%f %f %e\n", x, y, pt9[PRS->pnum(i,j)]);
		};
	};

	//flux x
	for(int i = 0; i <= 2 * gridN; i++)
	{
		for(int j = 0; j < gridN; j++)
		{
			x = i * PRS->h;
			y = PRS->h / 2.0 + j * PRS->h;
			fprintf(fd1, "%f %f %e\n", x, y, pt1[PRS->exnum(i,j)]);
			fprintf(fd9, "%f %f %e\n", x, y, pt9[PRS->exnum(i,j)]);
	
		};
	};

	//flux y
	for(int i = 0; i < 2 * gridN; i++)
	{
		for(int j = 0; j < gridN; j++)
		{
			x = PRS->h / 2.0 +i * PRS->h;
			y = j * PRS->h;

			fprintf(fd1, "%f %f %e\n", x, y, pt1[PRS->eynum(i,j)]);
			fprintf(fd9, "%f %f %e\n", x, y, pt9[PRS->eynum(i,j)]);
		};
	};

	fclose(fd1);
	fclose(fd9);

	delete[] pt1;
	delete[] pt9;
};

void printfluxmodule(solverparams * PRS, double * point, char * whereH, char * whereB)
{
	FILE * fdH, * fdB;
	double x, y;
	int gridN = PRS->n;
	double fx, fy, fxl, fxr, fyb, fyt;
	double hm;
	
	fdH = fopen(whereH, "w");
	fdB = fopen(whereB, "w");

	//values
	for(int i = 0; i < 2 * gridN; i++)
	{
		for(int j = 0; j < gridN; j++)
		{
			x = PRS->h / 2.0 + i * PRS->h;
			y = PRS->h / 2.0 + j * PRS->h;	

			fxl = point[PRS->exnum(i,j)];
			fxr = point[PRS->exnum(i + 1, j)];
			fyb = point[PRS->eynum(i, j)];
			fyt = point[PRS->eynum(i, (j+1) % gridN)];
			fx  = (fxl + fxr) / 2.0;
			fy  = (fyb + fyt) / 2.0;
			hm = sqrt(fx * fx + fy * fy);

			fprintf(fdH, "%f %f %e\n", x, y, hm);
			fprintf(fdB, "%f %f %e\n", x, y, PRS->cfunc(hm, x, y) * hm);

		};
		fprintf(fdH, "\n");
		fprintf(fdB, "\n");
	};
	fclose(fdH);
	fclose(fdB);
};

void printsimple(solverparams * PRS, double * point, char * wherev, char * wherex, char * wherey)
{
	FILE * fdv, * fdx, *fdy;
	double x, y;

	int gridN = PRS->n;

	fdv = fopen(wherev, "w");
	fdx = fopen(wherex, "w");
	fdy = fopen(wherey, "w");

	//values
	for(int i = 0; i < 2 * gridN; i++)
	{
		for(int j = 0; j < gridN; j++)
		{
			x = PRS->h / 2.0 + i * PRS->h;
			y = PRS->h / 2.0 + j * PRS->h;	

			fprintf(fdv, "%f %f %e\n", x, y, point[PRS->pnum(i,j)] * PRS->scaling);
		};
		fprintf(fdv, "\n");
	};

	//flux x
	for(int i = 0; i <= 2 * gridN; i++)
	{
		for(int j = 0; j < gridN; j++)
		{
			x = i * PRS->h;
			y = PRS->h / 2.0 + j * PRS->h;
		
			fprintf(fdx, "%f %f %e\n", x, y, point[PRS->exnum(i,j)]);
		};
		fprintf(fdx, "\n");
	};

	//flux y
	for(int i = 0; i < 2 * gridN; i++)
	{
		for(int j = 0; j < gridN; j++)
		{
			x = PRS->h / 2.0 +i * PRS->h;
			y = j * PRS->h;
			fprintf(fdy, "%f %f %e\n", x, y, point[PRS->eynum(i,j)]);
		};
		fprintf(fdy, "\n");
	};

};

void putsolution(solverparams * PRS, double * point)
{
	int gridN = PRS->n;

	double x,y;
	//values
	for(int i = 0; i < 2 * gridN; i++)
	{
		for(int j = 0; j < gridN; j++)
		{
			x = PRS->h / 2.0 + i * PRS->h;
			y = PRS->h / 2.0 + j * PRS->h;	
		
			point[PRS->pnum(i,j)] = PRS->solution(x,y);
		};
	};
	
	//flux x
	for(int i = 0; i <= 2 * gridN; i++)
	{
		for(int j = 0; j < gridN; j++)
		{
			x = i * PRS->h;
			y = PRS->h / 2.0 + j * PRS->h;
		
			point[PRS->exnum(i,j)] = PRS->fluxx(x,y);

		};
	};

	//flux y
	for(int i = 0; i < 2 * gridN; i++)
	{
		for(int j = 0; j < gridN; j++)
		{
			x = PRS->h / 2.0 +i * PRS->h;
			y = j * PRS->h;
			
			point[PRS->eynum(i,j)] = PRS->fluxy(x,y);
		};
	};
};

void printrs(solverparams * PRS, char * wherers)
{
	FILE * fdr;
	int gridN = PRS->n;
	double x,y, ys;
	fdr = fopen(wherers, "w");
	//rs
	for(int i = 0; i <= gridN; i++)
	{
		for(int j = 0; j < gridN; j++)
		{
			x = PRS->h / 2.0 + i * PRS->h;
			y = PRS->h / 2.0 + j * PRS->h;
				
			fprintf(fdr,"%f %f %e \n", x, y, PRS->rs(x,y));
		};
		fprintf(fdr,"\n");
	};
	for(int i = gridN + 1; i < 2 * gridN; i++)
	{
		for(int j = 0; j < gridN; j++)
		{
			x = PRS->h / 2.0 + i * PRS->h;
			y = PRS->h / 2.0 + j * PRS->h;
			ys = y + PRS->shiftnum * PRS->h;
			if (ys > 1) ys -=1;
			fprintf(fdr,"%f %f %e \n", x, y, PRS->rs(x,ys));
		};
		/*
		for(int j = PRS->n - PRS->shiftnum; j < PRS->n; j++)
		{
			x = PRS->h / 2.0 + i * PRS->h;
			y = PRS->h / 2.0 + j * PRS->h - 1;
			ys = y + PRS->shiftnum * PRS->h;
			fprintf(fdr,"%f %f %e \n", x, y, PRS->rs(x,ys));

		};
		
		for(int j = 0; j < PRS->n - PRS->shiftnum; j++)
		{
			x = PRS->h / 2.0 + i * PRS->h;
			y = PRS->h / 2.0 + j * PRS->h;
			ys = y + PRS->shiftnum * PRS->h;
			fprintf(fdr,"%f %f %e \n", x, y, PRS->rs(x,ys));

		};
		*/
		fprintf(fdr,"\n");
	};
	
	fclose(fdr);
};

void printreference(solverparams * PRS, double * point, double * refpoint, char * name)
{
	FILE * fdv, * fdx, * fdy;
	int gridN;
	double maxv, maxfx, maxfy;
	
	double * errorpoint;
	
	char namev[60];
	char namefx[60];
	char namefy[60];
	char nametv[60];
	char nametfx[60];
	char nametfy[60];
	char nameev[60];
	char nameefx[60];
	char nameefy[60];

	sprintf(namev, "%svalue.dat", name);
	sprintf(namefx, "%sfluxx.dat", name);
	sprintf(namefy, "%sfluxy.dat", name);

	sprintf(nametv, "%struevalue.dat", name);
	sprintf(nametfx, "%struefluxx.dat", name);
	sprintf(nametfy, "%struefluxy.dat", name);

	sprintf(nameev, "%svalueerror.dat", name);
	sprintf(nameefx, "%sfluxxerror.dat", name);
	sprintf(nameefy, "%sfluxyerror.dat", name);
	
	errorpoint = new double[PRS->varnum];
	gridN = PRS->n;
	
	maxv = 0; maxfx = 0; maxfy = 0;
	double tmp;

	for(int i = 0; i < PRS->varnum; i++)
	{
		errorpoint[i] = std::abs(point[i] - refpoint[i]);
	}; 

	for(int i = 0; i < 2 * gridN; i++)
	{
		for(int j = 0; j < gridN; j++)
		{
			tmp = std::abs(refpoint[PRS->pnum(i,j)]);	
			if (tmp > maxv) maxv = tmp;
		};
	};

	for(int i = 0; i < 2 * gridN + 1; i++)
	{
		for(int j = 0; j < gridN; j++)
		{
			tmp = std::abs(refpoint[PRS->exnum(i,j)]);	
			if (tmp > maxfx) maxfx = tmp;
		};
	};

	for(int i = 0; i < 2 * gridN; i++)
	{
		for(int j = 0; j < gridN; j++)
		{
			tmp = std::abs(refpoint[PRS->eynum(i,j)]);	
			if (tmp > maxfy) maxfy = tmp;
		};
	};
	
	for(int i = 0; i < 2 * gridN; i++)
	{
		for(int j = 0; j < gridN; j++)
		{
			errorpoint[PRS->pnum(i,j)] /= maxv * PRS->scaling;
		};
	};

	for(int i = 0; i < 2 * gridN + 1; i++)
	{
		for(int j = 0; j < gridN; j++)
		{
			errorpoint[PRS->exnum(i,j)] /= maxfx;
		};
	};

	for(int i = 0; i < 2 * gridN; i++)
	{
		for(int j = 0; j < gridN; j++)
		{
			errorpoint[PRS->eynum(i,j)] /= maxfy;
		};
	};

	printsimple(PRS, point, namev, namefx, namefy);
	printsimple(PRS, refpoint, nametv, nametfx, nametfy);
	printsimple(PRS, errorpoint, nameev, nameefx, nameefy);

	delete[] errorpoint;
};

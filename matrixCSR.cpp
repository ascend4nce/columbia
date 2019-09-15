#include "matrixCSR.h"
#include <mkl.h>
#include <mkl_rci.h>
#include <mkl_blas.h>
#include <mkl_spblas.h>
#include <mkl_service.h>
#include <cstdio>
#include <cmath>
#include "constants.h"

matrixCSR::matrixCSR()
{
	m = 0;
	n = 0;
	a = NULL; 
	ia = NULL;
	ja = NULL;
};

matrixCSR::~matrixCSR()
{
	if ((n > 0) && (m > 0))
	{
		delete[] a;
		delete[] ia;
		delete[] ja;
	};
};

//prints a sparse matrix in a console. For small matrices only.
void matrixCSR::print()
{
	int rowcounter, elementsinrow, nnzcounter, last, i, j;

	nnzcounter = 0;
	rowcounter = 0;

	printf("       ");
	for(i = 1; i <= n; i++)
	{
		printf(" <%3d> ", i);
	};
	printf("\n");

	while(rowcounter < m)
	{
		printf("  %3d|", rowcounter + 1);
		elementsinrow = ia[rowcounter + 1] - ia[rowcounter];
		last = 0;

		for(i = 0; i < elementsinrow; i++)
		{
			for(j = 0; j < ja[nnzcounter] - 1 - last; j++)
			{
				printf("       ");
			};

			printf(" % 02.3e", a[nnzcounter]);
			last = ja[nnzcounter++];
		};

		rowcounter++;
		printf("\n");
	};
};

matrixCSR::matrixCSR(int mset, int nset, double * aset, int * iaset, int * jaset)
{
	m = mset;
	n = nset;
	ia = new int[m+1];

	for(int i = 0; i < m + 1; i++)
	{
		ia[i] = iaset[i];
	}; 

	int nnz = iaset[m] - 1;

	ja = new int[nnz];
	a = new double[nnz];

	for(int i = 0; i < nnz; i++)
	{
		ja[i] = jaset[i];
		a[i] = aset[i];
	};
};

//Matrix-by-vector multiplication, right
void matrixCSR::matvec(double * vec, double * res)
{
	int  ctr;

	ctr = 0;
	for(int i = 0; i < m; i++)
	{
		res[i] = 0.0;

		for(int k = 0; k < ia[i + 1] - ia[i]; k++)
		{
			res[i] += vec[ja[ctr] - 1] * a[ctr];
			ctr++;
		};
	};
};

//Matrix-by-vector multiplication, left
void matrixCSR::leftmatvec(double * vec, double * res)
{
	int ctr = 0;

	for(int i = 0; i < n; i++)
	{
		res[i] = 0.0;
	};

	for(int i = 0; i < m; i++)
	{
		for(int k = 0; k < ia[i + 1] - ia[i]; k++)
		{
			res[ja[ctr] - 1] += vec[i] * a[ctr];
			ctr++;
		};
	};
};

//Driver for a GMRES implementation from Intel MKL, applied to a sparse matrix
void matrixCSR::solveprecgmres(double * x, double * rs, double lstolerance)
{
	int ione = 1;
	int izero = 0;
	double dzero = 0.0;
	double dminusone = -1.0;

	int ipar[128];
	double dpar[128];
	double * tmp;
	int iteration = 0;

	int sl_request;
	int ierr;
	int gmresit;

	double residual;
	
	double * trvec =  new double[n];
	double * trvec2 = new double[n];
	matrixCSR * precond;
	tmp = new double[(2 * gmrestart + 1) * n + gmrestart * (gmrestart + 9) + 1];
	//solving the system
	//preparations
	dcopy(&n, &dzero, &izero, x, &ione);
	dfgmres_init(&n, x, rs, &sl_request, ipar, dpar, tmp);

	ipar[4] = maxgmresit;
	ipar[10] = 1;

	ipar[7] = 1;
	ipar[8] = 1;
	ipar[9] = 0;
	dpar[0] = lstolerance;
	dpar[1] = 0.0;
	dpar[2] = dnrm2(&n, rs, &ione);
	dpar[3] = dpar[0] * dpar[2];
	ipar[2] = 1;
	ipar[3] = 0;

	dpar[30] = 1.0e-4;
	ipar[30] = 1;

	ipar[11] = 1;

	dpar[7] = 1.0e-14 * dnrm2(&n, rs, &ione);

	ipar[14] = gmrestart;

	int plen = sqrt(n / 6.0) * 2;
	int pcnnz = (2 * plen + 1) * n;
	double pctolerance = 1.0e-5;
	precond = new matrixCSR;
	precond->n = n;
	precond->m = n;
	precond->a = new double[pcnnz];
	precond->ia = new int[n + 1];
	precond->ja = new int[pcnnz];

	dcsrilut(&n, a, ia, ja, precond->a, precond->ia, precond->ja, &pctolerance, &plen, ipar, dpar, &ierr);


	dfgmres_check(&n, x, rs, &sl_request, ipar, dpar, tmp);
	sl_request = 100;
	
	int inum = 0;
	//running gmres
	while(sl_request > 0)
	{
		dfgmres(&n, x, rs, &sl_request, ipar, dpar, tmp);

		if (sl_request == 1)
		{
                        //printf("(%e %e)\n", dpar[3], dpar[4]);
			matvec(tmp + ipar[21] - 1, tmp + ipar[22] - 1);
		};
		
		if (sl_request == 3)
		{
			//printf("no pls");
			
			char cvar1, cvar, cvar2;
			cvar1 = 'L';
			cvar = 'N';
			cvar2 = 'U';
			mkl_dcsrtrsv(&cvar1, &cvar, &cvar2, &n, precond->a, precond->ia, precond->ja, tmp + ipar[21] - 1, trvec);

			cvar1 = 'U';
			cvar = 'N';
			cvar2 = 'N';
			mkl_dcsrtrsv(&cvar1, &cvar, &cvar2, &n, precond->a, precond->ia, precond->ja, trvec, tmp + ipar[22] - 1);
			
		};

		if((sl_request != 1) && (sl_request != 0) && (sl_request != 3))
		{
			printf("#Unexpected request %d\n", sl_request);
			break;
		};
	};
	if (sl_request != 0)
	{
		if (printlevel >= 1) printf("GMRES: Error %d encountered.\n", sl_request);
	};

	ipar[12] = 0;

	dfgmres_get(&n, x, rs, &sl_request, ipar, dpar, tmp, &gmresit);
	
	if (printlevel >= 2) printf("#GMRES it: %d# ", gmresit);

	delete precond;	
	delete[] trvec;
	delete[] trvec2;
	delete[] tmp;
}; 

//Driver for an accurate SVD-based linear system solution, applied to a sparse matrix; for the cases of ill-conditioned matrices
void matrixCSR::solvesvd(double * x, double * rs)
{
	double * matrix;
	double * tmpvec;

	double * U, * VT, * snum;

	if (m != n)
	{
		printf("Rectangular hessian.\n");
		return;
	};
	
	matrix = new double[n * n];
	tmpvec = new double[n];

	U = new double[n * n];
	VT = new double[n * n];
	snum = new double[n];

	int sz = n * n;
	int ione = 1;
	int izero = 0;
	double done = 1.0;
	double dzero = 0.0;
	char cA = 'A';
	char cT = 'T';
	int info;
	
	int lwork = 64 * n;
	double * work = new double[lwork];
	
	dcopy(&sz, &dzero, &izero, matrix, &ione);
	dcopy(&n, &dzero, &izero, tmpvec, &ione);
	

	for(int j = 0; j < n; j++)
	{
		if (j != 0)
		{
			tmpvec[j - 1] = 0.0;
		};
		tmpvec[j] = 1.0;

		matvec(tmpvec, matrix + j * n);
	};

	dgesvd(&cA, &cA, &n, &n, matrix, &n, snum, U, &n, VT, &n, work, &lwork, &info);
	if (info != 0) printf("dgesvd failed!\n");
/*
	for(int i = 0; i < n; i++)
	{
		printf("%e ", snum[i]);
	};
	printf("\n");
*/
	dgemv(&cT, &n, &n, &done, U, &n, rs, &ione, &dzero, tmpvec, &ione);
	for(int i = 0; i < n; i++)
	{
		if (snum[i] != 0)
		{
			tmpvec[i] /= snum[i];
		}
		else
		{
			if (printlevel >= 1) printf("SVD: zero singular value.\n");
			tmpvec[i] = 0.0;
		};
	};

	dgemv(&cT, &n, &n, &done, VT, &n, tmpvec, &ione, &dzero, x, &ione);
/*
	for(int i = 0; i < n; i++)
	{
		printf("%e ", x[i]);	
	};
	printf("\n");
	*/
	delete[] matrix;
	delete[] tmpvec;
	delete[] work;
	delete[] U;
	delete[] VT;
	delete[] snum;
};

//Driver for a CG-solution of regularized covariance matrix A^* A + \lambda I, applied to a sparse matrix
void matrixCSR::solvepreccovcg(double * x, double * rs, double lambda, double lstolerance)
{
	int ione = 1;
	int izero = 0;
	double dzero = 0.0;
	double done = 1.0;
	double dminuslambda = -lambda;

	int ipar[128];
	double dpar[128];

	double * temp = new double[n];
	double * tmp = new double[4 * n];
	double * trvec = new double[n];
	double * trvec2 = new double[n];

	int ierr;

	int sl_request = 0;

	dcopy(&n, &dzero, &izero, x, &ione);
	dcg_init(&n, x, rs, &sl_request, ipar, dpar, tmp);

	ipar[4] = maxcgit;
	ipar[10] = 1;

	ipar[7] = 1;
	ipar[8] = 1;
	ipar[9] = 0;
	dpar[1] = 0.0; //lstolerance * dnrm2(&n, rs, &ione);
	dpar[2] = 0.0;
	dpar[0] = lstolerance;
	dpar[3] = dpar[0] * dnrm2(&n, rs, &ione);
	ipar[2] = 1;
	ipar[3] = 0;

	//printf("!!%e!! \n\n", dnrm2(&n, rs, &ione));

	dcg_check(&n, x, rs, &sl_request, ipar, dpar, tmp);
	sl_request = 100;

	dpar[30] = 1.0e-4;
	ipar[30] = 1;

	matrixCSR * precond;
	
	int plen = sqrt(n / 6.0) * 2.0;
	int pcnnz = (2 * plen + 1) * n;
	double pctolerance = 1.0e-5;
	precond = new matrixCSR;
	precond->n = n;
	precond->m = n;
	precond->a = new double[pcnnz];
	precond->ia = new int[n + 1];
	precond->ja = new int[pcnnz];
	
	char cvar1, cvar, cvar2;
	
	dcsrilut(&n, a, ia, ja, precond->a, precond->ia, precond->ja, &pctolerance, &plen, ipar, dpar, &ierr);
	if (ierr != 0)
	{
		printf("Preconditioning failed.\n");
	};

	/*
	cvar1 = 'U';
	cvar = 'T';
	cvar2 = 'N';
	mkl_dcsrtrsv(&cvar1, &cvar, &cvar2, &n, precond->a, precond->ia, precond->ja, rs, trvec);

	cvar1 = 'L';
	cvar = 'T';
	cvar2 = 'U';
	mkl_dcsrtrsv(&cvar1, &cvar, &cvar2, &n, precond->a, precond->ia, precond->ja, trvec, rs);
	*/
	double vmax = 0.0;
/*
	for(int i = 0; i < n; i++)
	{
		double vlin = 0.0;
		for(int j = ia[i]; j < ia[i+1]; j++)
		{
			vlin += std::abs(a[j]);
		};
		if (vlin >= vmax) vmax = vlin;
	};
	
	printf("Regularization: %e %e\n.", vmax, lambda);	
*/	
	//running cg
	while(sl_request > 0)
	{
		dcg(&n, x, rs, &sl_request, ipar, dpar, tmp);


		if (sl_request == 1)
		{
			matvec(tmp, trvec);
			leftmatvec(trvec, tmp + n);
			daxpy(&n, &lambda, tmp, &ione, tmp + n, &ione);
		};

		if (sl_request == 3)
		{
			cvar1 = 'U';
			cvar = 'T';
			cvar2 = 'N';
			mkl_dcsrtrsv(&cvar1, &cvar, &cvar2, &n, precond->a, precond->ia, precond->ja, tmp + 2 * n, trvec);

			cvar1 = 'L';
			cvar = 'T';
			cvar2 = 'U';
			mkl_dcsrtrsv(&cvar1, &cvar, &cvar2, &n, precond->a, precond->ia, precond->ja, trvec, trvec2);

			cvar1 = 'L';
			cvar = 'N';
			cvar2 = 'U';
			mkl_dcsrtrsv(&cvar1, &cvar, &cvar2, &n, precond->a, precond->ia, precond->ja, trvec2, trvec);

			cvar1 = 'U';
			cvar = 'N';
			cvar2 = 'N';
			mkl_dcsrtrsv(&cvar1, &cvar, &cvar2, &n, precond->a, precond->ia, precond->ja, trvec, tmp + 3 * n);
			
		};

		if ((sl_request != 0) && (sl_request != 1) && (sl_request != 3))
		{
			if (printlevel >= 1) printf("CG: Error %d encountered.\n", sl_request);
		};
	};
	
	ipar[12] = 0;
	int cgit;
	dcg_get(&n, x, rs, &sl_request, ipar, dpar, tmp, &cgit);

	if (printlevel >= 2) printf("(#CG it: %d#)", cgit);
	delete precond;
	delete[] temp;
	delete[] tmp;
	delete[] trvec;
	delete[] trvec2;
};

//Driver for accurate SVD-based solution of a linear system with a regularized covariance matrix A^* A + \lambda I, for ill-conditioned cases
void matrixCSR::solvecovsvd(double * x, double * rs, double lambda)
{
	double * cov;
	double * matrix;
	double * tmpvec;

	double * U, * VT, * snum;

	if (m != n)
	{
		printf("Rectangular hessian.\n");
		return;
	};
	
	matrix = new double[n * n];
	cov = new double[n * n];	

	tmpvec = new double[n];

	U = new double[n * n];
	VT = new double[n * n];
	snum = new double[n];

	int sz = n * n;
	int ione = 1;
	int izero = 0;
	double done = 1.0;
	double dzero = 0.0;
	char cA = 'A';
	char cN = 'N';
	char cT = 'T';
	int info;
	
	int lwork = 64 * n;
	double * work = new double[lwork];
	
	dcopy(&sz, &dzero, &izero, matrix, &ione);
	dcopy(&n, &dzero, &izero, tmpvec, &ione);
	dcopy(&sz, &dzero, &izero, cov, &ione);	

	for(int j = 0; j < n; j++)
	{
		cov[j * (n + 1)] = lambda;
		if (j != 0)
		{
			tmpvec[j - 1] = 0.0;
		};
		tmpvec[j] = 1.0;

		matvec(tmpvec, matrix + j * n);
	};

	dgemm(&cT, &cN, &n, &n, &n, &done, matrix, &n, matrix, &n, &done, cov, &n);

	dgesvd(&cA, &cA, &n, &n, cov, &n, snum, U, &n, VT, &n, work, &lwork, &info);
	
	//printf("\n(%e )", snum[n-1]);
	//printf("\n");

	if (info != 0) printf("dgesvd failed!\n");

	dgemv(&cT, &n, &n, &done, U, &n, rs, &ione, &dzero, tmpvec, &ione);
	for(int i = 0; i < n; i++)
	{
		if (snum[i] != 0)
		{
			tmpvec[i] /= snum[i];
		}
		else
		{
			if (printlevel>=1) printf("SVD: zero singular value.\n");
			tmpvec[i] = 0.0;
		};
	};

	dgemv(&cT, &n, &n, &done, VT, &n, tmpvec, &ione, &dzero, x, &ione);

	delete[] matrix;
	delete[] tmpvec;
	delete[] work;
	delete[] U;
	delete[] VT;
	delete[] snum;
	delete[] cov;
};



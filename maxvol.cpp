//This code is an implementation of a greedy Maximum-Volume-Submatrix algorithm to a matrix with unitary columns.

#include <cstdlib>
#include <cmath>
#include <mkl.h>
#include <mkl_blas.h>
#include <cstdio>

int maxvol2(int r, int n, const double *factor, int lda, int *J)
{
    const int size = n*r;
    int i, j;
    const int ione = 1;
    int *tmp_ipiv;
	int info;
    double *tmp_matrix, *tmp_row, *tmp_column, *pointer1, *pointer2;
    double alpha = 1.0, maxf;
    tmp_matrix = (double *)malloc(r*r*sizeof(double));
    tmp_ipiv = (int *)malloc(r*sizeof(int));

	double * coef = new double[r * n];

	for(int j = 0; j < r; j++)
	{
		for(int i = 0; i < n; i++)
		{
			coef[j + i * r] = factor[i + j * lda];
		};
	};

	double eps = 1.0e-8;
    //dcopy_(&size, factor, &ione, coef, &ione);
    for (i = 0; i < r; i++)
    {
        dcopy(&r, coef + r * J[i], &ione, tmp_matrix + i*r, &ione);
    }
    dgesv(&r, &n, tmp_matrix, &r, tmp_ipiv, coef, &r, &info);
	if(info != 0) printf("Maxvol: problem (%d).\n", info);
    free(tmp_ipiv);
    free(tmp_matrix);
    tmp_row = (double *) malloc(n*sizeof(double));
    tmp_column = (double *) malloc(r*sizeof(double));
    maxf = eps + 2;
    while (std::abs(maxf) > eps + 1.0)
    {
        i = idamax(&size, coef, &ione) - 1;
        maxf = coef[i];
        j = i/r;
        i -= j*r;
        if (std::abs(maxf) > eps + 1.0)
        {
            dcopy(&n, coef+i, &r, tmp_row, &ione);
            dcopy(&r, coef+j*r, &ione, tmp_column, &ione);
            tmp_column[i] -= 1.0;
            int k;
            for (k = 0; k < n; k++)
            {
                if (J[k] == j)
                {
                    break;
                }
            }
            J[k] = J[i];
            J[i] = j;
		
            alpha = -1.0/maxf;
            dger(&r, &n, &alpha, tmp_column, &ione, tmp_row, &ione, coef, &r);
        }
    }
    free(tmp_row);
    free(tmp_column);

	delete[] coef;
    return 0;
}

int maxvol1(int r, int n, double *a,int lda, int * J)
{
    for (int i = 0; i < n; i++)
    {
        J[i] = i;
    }
    maxvol2(r, n, a, lda, J);
}



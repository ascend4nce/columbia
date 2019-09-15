#include <mkl.h>
#include <mkl_blas.h>
#include <cmath>
#include <cstdio>
#include "matrixCSR.h"
#include "nonlinearities.h"
#include "fullscheme.h"

int modL(int num, int loop)
{
	if ((num >= 0) && (num < loop)) return num;

	if (num < 0) return num + loop;

	if (num >= loop) return num - loop;
};

solverparams::solverparams(int nset, double alphaset, double betaset, double shiftxs, double shiftys, double scalingset)
{
	alpha = alphaset;
	beta = betaset;
	scaling = scalingset;

	n = nset;
	h = 1.0 / (double) nset;

	varnum = 6 * n * n + n;

	nlnum = 4 * n * n + n;

	halfnum = std::floor((shiftxs / h) + 0.00001);
	shiftnum = std::floor((shiftys / h) + 0.00001);
	
	shifty = shiftnum * h;
	shiftx = halfnum * h;
};	


double solverparams::cfunc(double s, double x, double y)
{
	//return test12NL(s, x, y, alpha, beta, shiftx, shifty);
	return SimpleNL(s, x, y, alpha, beta, shiftx, shifty);
	//return MagnetNL(s, x, y, alpha, beta, shiftx, shifty);
	//return EngineNL(s, x, y, shiftx, shifty);
};

double solverparams::cfuncderiv(double s, double x, double y)
{
	//return test12Deriv(s,x,y, alpha, beta,shiftx, shifty);
	return SimpleDeriv(s,x,y, alpha, beta,shiftx, shifty);
	//return MagnetDeriv(s,x,y, alpha, beta,shiftx, shifty);
	//return EngineDeriv(s,x,y, shiftx, shifty);
};

double solverparams::rs(double x, double y)
{
	//return test12Rs(x,y, alpha, beta, shiftx, shifty);
	return SimpleRs(x,y, alpha, beta, shiftx, shifty);
	//return EngineRs(x,y, beta, shiftx, shifty);
	//return MagnetRs(x,y, alpha, beta, shiftx, shifty);
};

double solverparams::solution(double x, double y)
{
	return 0;
};

double solverparams::fluxx(double x, double y)
{
	return 0;
};

double solverparams::fluxy(double x, double y)
{
	return 0;
};

//auxillary functions defining the memory positions of values, x and y fluxes in the vector of unknowns
int solverparams::pnum(int i, int j)
{
	return 4 * n * n + n + i * n + j;
};

int solverparams::exnum(int i, int j)
{
	return i * n + j;
};

int solverparams::eynum(int i, int j)
{
	return 2 * n * n + n + n * i + j;
};

int solverparams::shiftindex(int i)
{
	int xc, yc;
	yc = i % n;
	xc = i / n;

	if (((xc >= halfnum + 1) && (xc < 2 * n + 1)) || ((xc >= 2 * n + halfnum + 1) && (xc < 4 * n + 1)) || ((xc >= 4 * n + halfnum + 1) && (xc < 6 * n + 1)))
	{
		yc = modL(yc - shiftnum, n);
	};

	return yc + xc * n;
};

int solverparams::unrollindex(int i)
{
	int xc, yc;
	yc = i % n;
	xc = i / n;

	if (((xc >= halfnum + 1) && (xc < 2 * n + 1)) || ((xc >= 2 * n  + halfnum + 1) && (xc < 4 * n + 1)) || ((xc >= 4 * n + halfnum + 1) && (xc < 6 * n + 1)))
	{
		yc = modL(yc + shiftnum, n);
	};

	return yc + xc * n;
};


//auxilary function that aligns shifted solution vectors
void solverparams::vecshiftrotate(double * values)
{	
	double * tmp;
	tmp = new double[varnum];

	for(int i = 0; i < varnum; i++)
	{
		tmp[shiftindex(i)] = values[i];
	};

	for(int i = 0; i < varnum; i++)
	{
		values[i] = tmp[i];
	};

	delete[] tmp;
};

//auxilary function that aligns shifted solution vectors
void solverparams::vecshiftunroll(double * values)
{	
	double * tmp;
	tmp = new double[varnum];

	for(int i = 0; i < varnum; i++)
	{
		tmp[unrollindex(i)] = values[i];
	};

	for(int i = 0; i < varnum; i++)
	{
		values[i] = tmp[i];
	};

	delete[] tmp;
};



double oneresidlin(double * pt, void * params, int num)
{
	solverparams * prs;
	prs = (solverparams *) params;

	int i, j;
	double res = 0.0;
	
	if (num < (2 * prs->n + 1) * prs->n)
	{
		i = num / prs->n;
		j = num - i * prs->n;

		//nl equation, x
		if (i != 0)
		{
			res -=  -pt[prs->pnum(i - 1, j)] * prs->scaling / prs->h;
		};
		if (i != 2 * prs->n)
		{
			res -= pt[prs->pnum(i, j)] * prs->scaling / prs->h;
		};
	}
	else
	{
		num -= (2 * prs->n + 1) * prs->n;
		if (num < prs->n * 2 * prs->n)
		{
			i = num / prs->n;
			j = num - i * prs->n;
			
			//nl equation, y
			if (j == 0)
                       	{
				res -= (pt[prs->pnum(i, 0)] - pt[prs->pnum(i, prs->n - 1)]) * prs->scaling / prs->h;
                       	}
			else
			{
				res -= (pt[prs->pnum(i, j)] - pt[prs->pnum(i, j - 1)]) * prs->scaling / prs->h;
			};                  
			
		}
		else
		{
			num -= 2 * prs->n * prs->n;
			
			i = num / prs->n;
			j = num - i * prs->n;

			//simple equation
			res += (pt[prs->exnum(i + 1, j)] - pt[prs->exnum(i,j)]) * prs->scaling / prs->h;
			res += (pt[prs->eynum(i, modL(j + 1, prs->n))] - pt[prs->eynum(i,j)]) * prs->scaling / prs->h;
		};
	};

	return res;
};

double oneresidnl(double * pt, void * params, int num)
{
	solverparams * prs;
	prs = (solverparams *) params;

	int i, j;
	double res = 0.0;
	
	double x, y;
	
	if (num < (2 * prs->n + 1) * prs->n)
	{
		i = num / prs->n;
		j = num - i * prs->n;

		//nl equation, x
		double tmp[4];
		double nlt, nlb, nrt, nrb, nlc, nrc;
				
		// nonlinear part
		//left
		x = prs->h / 2.0 + prs->h * (i - 1);
		y = prs->h / 2.0 + prs->h * j;
                if (i != 0)
            	{
                        tmp[0] = pt[prs->exnum(i, j)];
                        tmp[3] = pt[prs->exnum(i - 1, j)];

                        tmp[1] = pt[prs->eynum(i - 1, modL(j + 1, prs->n))];
                        tmp[2] = pt[prs->eynum(i - 1, j)];

                        nlt = sqrt(tmp[1] * tmp[1] + 0.25 * (tmp[3] + tmp[0]) * (tmp[3] + tmp[0]));
                        nlb = sqrt(tmp[2] * tmp[2] + 0.25 * (tmp[3] + tmp[0]) * (tmp[3] + tmp[0]));
                        nlc = sqrt(tmp[0] * tmp[0] + 0.25 * (tmp[1] + tmp[2]) * (tmp[1] + tmp[2]));

			res += (tmp[0] + tmp[3]) * (prs->cfunc(nlt, x, y) + prs->cfunc(nlb, x, y)) / 16.0;
                        res += tmp[0] * prs->cfunc(nlc, x, y) / 4.0;
		};
		x += prs->h;
                //right
                if ((i != 2 * prs->n))
                {
                        tmp[0] = pt[prs->exnum(i, j)];
                        tmp[3] = pt[prs->exnum(i + 1, j)];

                        tmp[1] = pt[prs->eynum(i, modL(j + 1, prs->n))];
                        tmp[2] = pt[prs->eynum(i, j)];

                        nrt = sqrt(tmp[1] * tmp[1] + 0.25 * (tmp[3] + tmp[0]) * (tmp[3] + tmp[0]));
                        nrb = sqrt(tmp[2] * tmp[2] + 0.25 * (tmp[3] + tmp[0]) * (tmp[3] + tmp[0]));
                        nrc = sqrt(tmp[0] * tmp[0] + 0.25 * (tmp[1] + tmp[2]) * (tmp[1] + tmp[2]));
                                
                        res += (tmp[0] + tmp[3]) * (prs->cfunc(nrt,x,y) + prs->cfunc(nrb,x,y)) / 16.0;
                        res += tmp[0] * prs->cfunc(nrc,x,y) / 4.0;
                };
	}
	else
	{
		num -= (2 * prs->n + 1) * prs->n;
		if (num < (prs->n) * 2 * prs->n)
		{
			i = num / (prs->n);
			j = num - i * (prs->n);
			
			//nl equation, y
			double tmp[4];
			double ntc, ntl, ntr, nbc, nbl, nbr;

			//nonlinear part
			//top
			x = prs->h / 2.0 + prs->h * i;
			y = prs->h / 2.0 + prs->h * j;

			tmp[0] = pt[prs->eynum(i, j)];
			tmp[3] = pt[prs->eynum(i, modL(j + 1, prs->n))];

			tmp[1] = pt[prs->exnum(i, j)];
			tmp[2] = pt[prs->exnum(i + 1, j)];

			ntl = sqrt(tmp[1] * tmp[1] + 0.25 * (tmp[3] + tmp[0]) * (tmp[3] + tmp[0]));
			ntr = sqrt(tmp[2] * tmp[2] + 0.25 * (tmp[3] + tmp[0]) * (tmp[3] + tmp[0]));
			ntc = sqrt(tmp[0] * tmp[0] + 0.25 * (tmp[1] + tmp[2]) * (tmp[1] + tmp[2]));

			res += (tmp[0] + tmp[3]) * (prs->cfunc(ntl,x,y) + prs->cfunc(ntr,x,y)) / 16.0;
			res += tmp[0] * prs->cfunc(ntc,x,y) / 4.0;

			//bottom
			y -= prs->h;
			if (y < 0) 
			{
				y += 1.0;
			};

			tmp[0] = pt[prs->eynum(i, j)];
			tmp[3] = pt[prs->eynum(i, modL(j - 1, prs->n))];

			tmp[1] = pt[prs->exnum(i, modL(j - 1, prs->n))];
			tmp[2] = pt[prs->exnum(i + 1, modL(j - 1, prs->n))];

			nbl = sqrt(tmp[1] * tmp[1] + 0.25 * (tmp[3] + tmp[0]) * (tmp[3] + tmp[0]));
			nbr = sqrt(tmp[2] * tmp[2] + 0.25 * (tmp[3] + tmp[0]) * (tmp[3] + tmp[0]));
			nbc = sqrt(tmp[0] * tmp[0] + 0.25 * (tmp[1] + tmp[2]) * (tmp[1] + tmp[2]));

			res += (tmp[0] + tmp[3]) * (prs->cfunc(nbl,x,y) + prs->cfunc(nbr,x,y)) / 16.0;
			res += tmp[0] * prs->cfunc(nbc,x,y) / 4.0;
		}
		else
		{
			num -= 2 * prs->n * (prs->n);
			
			i = num / prs->n;
			j = num - i * prs->n;

			//simple equation
			res -= prs->rs(i * prs->h + prs->h / 2.0, j * prs->h + prs->h / 2.0) * prs->scaling;
		};
	};

	return res;
};

void fullsolverresidual(double * pt, double * ev, void * params)
{
	solverparams * prs;
	prs = (solverparams *) params;

	double res;
	for(int i = 0; i < prs->varnum; i++)
	{
		res = 0.0;
		res += oneresidnl(pt, params, i);
		res += oneresidlin(pt, params, i);
		ev[i] = res;
	};
};

matrixCSR * fullsolverderivnum(double * pt, void * params)
{
	solverparams * prs;
	prs = (solverparams *) params;

        double * tempup, * tempdown, * temppoint;
	int ione = 1;
        int izero = 0;
	double done = 1.0;
        double dzero = 0.0;
	double dminusone = -1.0;

	int vn = prs->varnum;
	double h = 1e-7;
        double hinv = 1.0 / h;

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
	        temppoint[j] += h / 2.0;
                fullsolverresidual(temppoint, tempup, params);
                temppoint[j] -= h;
                fullsolverresidual(temppoint, tempdown, params);
                temppoint[j] += h / 2.0;

                daxpy(&vn, &dminusone, tempdown, &ione, tempup, &ione);
                dscal(&vn, &hinv, tempup, &ione);

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


matrixCSR * fullsolverderivCSR(double * pt, void * params)
{
	solverparams * prs;
	prs = (solverparams *) params;

	matrixCSR * result;

	result    = new matrixCSR();

	int vn = prs->varnum;
	int nnz = 44 * prs->n * prs->n + 1 * prs->n;
	
	result->n = vn;
	result->m = vn;
	result->a = new double[nnz];
        result->ia = new int[vn + 1];
        result->ja = new int[nnz];


	int ctr = 0;
	int eqnum = 0;
	result->ia[0] = ctr;
	
	double x,y;

	//nl equation, x
	for(int i = 0; i <= 2 * prs->n; i++)
	{
		for(int j = 0; j < prs-> n; j++)
		{
			double tmp[8];
			double nlt, nlb, nlc, nrt, nrb, nrc;
			double res;
			if (i != 0)
			{
				tmp[0] = pt[prs->exnum(i, j)];
				tmp[3] = pt[prs->exnum(i - 1, j)];
		
				tmp[1] = pt[prs->eynum(i - 1, j)];
				tmp[2] = pt[prs->eynum(i - 1, modL(j + 1, prs->n))];
			};
			
			if (i != 2 * prs->n)
			{
				tmp[4] = pt[prs->exnum(i, j)];
				tmp[7] = pt[prs->exnum(i + 1, j)];
			
				tmp[5] = pt[prs->eynum(i, j)];
				tmp[6] = pt[prs->eynum(i, modL(j+1, prs->n))];
			};
			
			if (i != 0)
			{
				nlc = sqrt(tmp[0] * tmp[0] + 0.25 * (tmp[1] + tmp[2]) * (tmp[1] + tmp[2]));
				nlb = sqrt(tmp[1] * tmp[1] + 0.25 * (tmp[0] + tmp[3]) * (tmp[0] + tmp[3]));
				nlt = sqrt(tmp[2] * tmp[2] + 0.25 * (tmp[0] + tmp[3]) * (tmp[0] + tmp[3]));
			};

			if (i != 2 * prs->n)
			{
				nrc = sqrt(tmp[4] * tmp[4] + 0.25 * (tmp[5] + tmp[6]) * (tmp[5] + tmp[6]));
				nrb = sqrt(tmp[5] * tmp[5] + 0.25 * (tmp[4] + tmp[7]) * (tmp[4] + tmp[7]));
				nrt = sqrt(tmp[6] * tmp[6] + 0.25 * (tmp[4] + tmp[7]) * (tmp[4] + tmp[7]));
			};
				
			//3
			x = prs->h/2.0 + (i-1) * prs->h;
			y = prs->h/2.0 + j * prs->h;

			if (i != 0)
			{
				res = 0;
				result->ja[ctr] = prs->exnum(i-1, j);
					
				res+= (prs->cfunc(nlt,x,y) + prs->cfunc(nlb,x,y)) / 16.0;
					
				res+= (0.25 * prs->cfuncderiv(nlt,x,y) * (tmp[0] + tmp[3]) / nlt) * (tmp[0] + tmp[3]) / 16.0;
				res+= (0.25 * prs->cfuncderiv(nlb,x,y) * (tmp[0] + tmp[3]) / nlb) * (tmp[0] + tmp[3]) / 16.0;
				result->a[ctr] = res;
				ctr++;
			};

			//4
			res = 0;
			result->ja[ctr] = prs->exnum(i, j);
			//left
			x = prs->h/2.0 + (i-1) * prs->h;
			y = prs->h/2.0 + j * prs->h;

			if (i != 0)
			{
				res += prs->cfunc(nlt,x,y) / 16.0 + prs->cfunc(nlb,x,y) / 16.0 + prs->cfunc(nlc,x,y) / 4.0;
					
				res += (prs->cfuncderiv(nlc,x,y) * tmp[0] / nlc) * tmp[0] / 4.0;
				res += (0.25 * prs->cfuncderiv(nlb,x,y) * (tmp[0] + tmp[3]) / nlb) * (tmp[0] + tmp[3]) / 16.0;
				res += (0.25 * prs->cfuncderiv(nlt,x,y) * (tmp[0] + tmp[3]) / nlt) * (tmp[0] + tmp[3]) / 16.0;
			};
			
			//right
			x = prs->h/2.0 + i * prs->h;
			y = prs->h/2.0 + j * prs->h;

			if (i != 2 * prs->n)
			{
				res += prs->cfunc(nrt,x,y) / 16.0 + prs->cfunc(nrb,x,y) / 16.0 + prs->cfunc(nrc,x,y) / 4.0;
					
				res += (prs->cfuncderiv(nrc,x,y) * tmp[4] / nrc) * tmp[4] / 4.0;
				res += (0.25 * prs->cfuncderiv(nrb,x,y) * (tmp[4] + tmp[7]) / nrb) * (tmp[4] + tmp[7]) / 16.0;
				res += (0.25 * prs->cfuncderiv(nrt,x,y) * (tmp[4] + tmp[7]) / nrt) * (tmp[4] + tmp[7]) / 16.0;
			};
			result->a[ctr] = res;
			ctr++;
	
			//5
			x = prs->h/2.0 + i * prs->h;
			y = prs->h/2.0 + j * prs->h;
			if (i != 2 * prs->n)
			{
				res = 0;
				result->ja[ctr] = prs->exnum(i + 1, j);

				res+= (prs->cfunc(nrt,x,y) + prs->cfunc(nrb,x,y)) / 16.0;
				
				res+= (0.25 * prs->cfuncderiv(nrt,x,y) * (tmp[4] + tmp[7]) / nrt) * (tmp[4] + tmp[7]) / 16.0;
				res+= (0.25 * prs->cfuncderiv(nrb,x,y) * (tmp[4] + tmp[7]) / nrb) * (tmp[4] + tmp[7]) / 16.0;
				result->a[ctr] = res;
				ctr++;
			};

			//6
			x = prs->h/2.0 + (i - 1) * prs->h;
			y = prs->h/2.0 + j * prs->h;
			if (i != 0)
			{
				if (j != prs->n - 1)
				{
					res = 0;
					result->ja[ctr] = prs->eynum(i - 1, j);

					res += (prs->cfuncderiv(nlb,x,y) * tmp[1] / nlb) * (tmp[0] + tmp[3]) / 16.0;
					res += (prs->cfuncderiv(nlc,x,y) * 0.25 * (tmp[1] + tmp[2]) / nlc) * tmp[0] / 4.0;

					result->a[ctr] = res;
					ctr++;
				};
			};

			//7
			if (i != 0)
			{
				res = 0;
				result->ja[ctr] = prs->eynum(i - 1, modL(j + 1, prs->n));

				res += (prs->cfuncderiv(nlt,x,y) * tmp[2] / nlt) * (tmp[0] + tmp[3]) / 16.0;
				res += (prs->cfuncderiv(nlc,x,y) * 0.25 * (tmp[1] + tmp[2]) / nlc) * tmp[0] / 4.0;
			
				result->a[ctr] = res;
				ctr++;
			};

			//6, loop
			if (i != 0)
			{
				if (j == prs->n - 1)
				{
					res = 0;
					result->ja[ctr] = prs->eynum(i - 1, j);

					res += (prs->cfuncderiv(nlb,x,y) * tmp[1] / nlb) * (tmp[0] + tmp[3]) / 16.0;
					res += (prs->cfuncderiv(nlc,x,y) * 0.25 * (tmp[1] + tmp[2]) / nlc) * tmp[0] / 4.0;

					result->a[ctr] = res;
					ctr++;
				};
			};

			//8
			x = prs->h/2.0 + i * prs->h;
			y = prs->h/2.0 + j * prs->h;

			if (i != 2 * prs->n)
			{
				if (j != prs->n - 1)
				{
					res = 0;
					result->ja[ctr] = prs->eynum(i, j);

					res += (prs->cfuncderiv(nrb,x,y) * tmp[5] / nrb) * (tmp[4] + tmp[7]) / 16.0;
					res += (prs->cfuncderiv(nrc,x,y) * 0.25 * (tmp[5] + tmp[6]) / nrc) * tmp[4] / 4.0;

					result->a[ctr] = res;
					ctr++;
				};
			};

			//9
			if (i != 2 * prs->n)
			{
				res = 0;
				result->ja[ctr] = prs->eynum(i, modL(j + 1, prs->n));

				res += (prs->cfuncderiv(nrt,x,y) * tmp[6] / nrt) * (tmp[4] + tmp[7]) / 16.0;
				res += (prs->cfuncderiv(nrc,x,y) * 0.25 * (tmp[5] + tmp[6]) / nrc) * tmp[4] / 4.0;
				
				result->a[ctr] = res;
				ctr++;
			};

			//8, loop
			if (i != 2 * prs->n)
			{
				if (j == prs->n - 1)
				{
					res = 0;
					result->ja[ctr] = prs->eynum(i, j);

					res += (prs->cfuncderiv(nrb,x,y) * tmp[5] / nrb) * (tmp[4] + tmp[7]) / 16.0;
					res += (prs->cfuncderiv(nrc,x,y) * 0.25 * (tmp[5] + tmp[6]) / nrc) * tmp[4] / 4.0;

					result->a[ctr] = res;
					ctr++;
				};
			};

			//gradientX
			//1
			if (i != 0)
			{
				result->ja[ctr] = prs->pnum(i - 1, j);
				result->a[ctr] = + 1.0 * prs->scaling / prs->h;
				ctr++;
			};
			//2
			if (i != 2 * prs->n)
			{
				result->ja[ctr] = prs->pnum(i, j);
				result->a[ctr] = - 1.0 * prs->scaling / prs->h;
				ctr++;
			};
			
			result->ia[eqnum + 1] = ctr;
			eqnum++;
		};
	};

	//nl equation, y
	for(int i = 0; i < 2 * prs->n; i++)
	{
		for(int j = 0; j < prs->n; j++)
		{
			double tmp[8];
			double ntc, ntl, ntr, nbc, nbl, nbr;
			double res;
			
			tmp[0] = pt[prs->eynum(i, j)];
			tmp[3] = pt[prs->eynum(i, modL(j - 1, prs->n))];
		
			tmp[1] = pt[prs->exnum(i, modL(j - 1, prs->n))];
			tmp[2] = pt[prs->exnum(i + 1, modL(j - 1, prs->n))];
				
			tmp[4] = pt[prs->eynum(i, j)];
			tmp[7] = pt[prs->eynum(i, modL(j + 1, prs->n))];
			
			tmp[5] = pt[prs->exnum(i, j)];
			tmp[6] = pt[prs->exnum(i + 1, j)];

			nbc = sqrt(tmp[0] * tmp[0] + 0.25 * (tmp[1] + tmp[2]) * (tmp[1] + tmp[2]));
			nbl = sqrt(tmp[1] * tmp[1] + 0.25 * (tmp[0] + tmp[3]) * (tmp[0] + tmp[3]));
			nbr = sqrt(tmp[2] * tmp[2] + 0.25 * (tmp[0] + tmp[3]) * (tmp[0] + tmp[3]));

			ntc = sqrt(tmp[4] * tmp[4] + 0.25 * (tmp[5] + tmp[6]) * (tmp[5] + tmp[6]));
			ntl = sqrt(tmp[5] * tmp[5] + 0.25 * (tmp[4] + tmp[7]) * (tmp[4] + tmp[7]));
			ntr = sqrt(tmp[6] * tmp[6] + 0.25 * (tmp[4] + tmp[7]) * (tmp[4] + tmp[7]));

			//3
			if (j != 0)
			{
				res = 0;

				x = prs->h/2.0 + i * prs->h;
				y = prs->h/2.0 + (j - 1) * prs->h;

				if (y < 0) y+= 1.0;

				result->ja[ctr] = prs->exnum(i, modL(j - 1, prs->n));
			
				res += (prs->cfuncderiv(nbl,x,y) * tmp[1] / nbl) * (tmp[0] + tmp[3]) / 16.0;
				res += (prs->cfuncderiv(nbc,x,y) * 0.25 * (tmp[1] + tmp[2]) / nbc) * tmp[0] / 4.0;
			
				result->a[ctr] = res;
				ctr++;
			};

			//4
			x = prs->h/2.0 + i * prs->h;
			y = prs->h/2.0 + j * prs->h;

			res = 0;
			result->ja[ctr] = prs->exnum(i, j);
				
			res += (prs->cfuncderiv(ntl,x,y) * tmp[5] / ntl) * (tmp[4] + tmp[7]) / 16.0;
			res += (prs->cfuncderiv(ntc,x,y) * 0.25 * (tmp[5] + tmp[6]) / ntc) * tmp[4] / 4.0;

			result->a[ctr] = res;
			ctr++;

			//3.5
			if (j == 0)
			{
				x = prs->h/2.0 + i * prs->h;
				y = prs->h/2.0 + (j - 1) * prs->h;

				if (y < 0) y+= 1.0;
				res = 0;

				result->ja[ctr] = prs->exnum(i, modL(j-1, prs->n));
			
				res += (prs->cfuncderiv(nbl,x,y) * tmp[1] / nbl) * (tmp[0] + tmp[3]) / 16.0;
				res += (prs->cfuncderiv(nbc,x,y) * 0.25 * (tmp[1] + tmp[2]) / nbc) * tmp[0] / 4.0;
				
				result->a[ctr] = res;
				ctr++;
			};


			//5
			if (j != 0)
			{
				x = prs->h/2.0 + i * prs->h;
				y = prs->h/2.0 + (j-1) * prs->h;
				if (y < 0) y+=1.0;

				res =  0;
				result->ja[ctr] = prs->exnum(i+1, modL(j - 1, prs->n));
			
				res += (prs->cfuncderiv(nbr,x,y) * tmp[2] / nbr) * (tmp[0] + tmp[3]) / 16.0;
				res += (prs->cfuncderiv(nbc,x,y) * 0.25 * (tmp[1] + tmp[2]) / nbc) * tmp[0] / 4.0;
				
				result->a[ctr] = res;
				ctr++;
			};

			//6
			res = 0;
			result->ja[ctr] = prs->exnum(i + 1, j);
			x = prs->h/2.0 + i * prs->h;
			y = prs->h/2.0 + j * prs->h;

			res += (prs->cfuncderiv(ntr,x,y) * tmp[6] / ntr) * (tmp[4] + tmp[7]) / 16.0;
			res += (prs->cfuncderiv(ntc,x,y) * 0.25 * (tmp[5] + tmp[6]) / ntc) * tmp[4] / 4.0;
	

			result->a[ctr] = res;
			ctr++;

			//5.5
			if (j == 0)
			{
				res =  0;
	
				x = prs->h/2.0 + i * prs->h;
				y = prs->h/2.0 + (j-1) * prs->h;

				if (y < 0) y+=1.0;
				result->ja[ctr] = prs->exnum(i+1, modL(j - 1, prs->n));
			
				res += (prs->cfuncderiv(nbr,x,y) * tmp[2] / nbr) * (tmp[0] + tmp[3]) / 16.0;
				res += (prs->cfuncderiv(nbc,x,y) * 0.25 * (tmp[1] + tmp[2]) / nbc) * tmp[0] / 4.0;
				
				result->a[ctr] = res;
				ctr++;
			};

			//9, loop
			if (j == prs->n - 1)
			{
				res = 0;
				result->ja[ctr] = prs->eynum(i, modL(j + 1, prs->n));
				x = prs->h/2.0 + i * prs->h;
				y = prs->h/2.0 + j * prs->h;
				if (y > 1.0) y -= 1.0;

				res+= (prs->cfunc(ntl,x,y) + prs->cfunc(ntr,x,y)) / 16.0;

				res+= (0.25 * prs->cfuncderiv(ntl,x,y) * (tmp[4] + tmp[7]) / ntl) * (tmp[4] + tmp[7]) / 16.0;
				res+= (0.25 * prs->cfuncderiv(ntr,x,y) * (tmp[4] + tmp[7]) / ntr) * (tmp[4] + tmp[7]) / 16.0;

				result->a[ctr] = res;
				ctr++;
			};

			//7			
			if (j != 0)
			{
				x = prs->h/2.0 + i * prs->h;
				y = prs->h/2.0 + (j - 1) * prs->h;
				if (y < 0) y += 1.0;

				res = 0;
				result->ja[ctr] = prs->eynum(i, modL(j - 1, prs->n));

				res+= (prs->cfunc(nbl,x,y) + prs->cfunc(nbr,x,y)) / 16.0;

				res+= (0.25 * prs->cfuncderiv(nbl,x,y) * (tmp[0] + tmp[3]) / nbl) * (tmp[0] + tmp[3]) / 16.0;
				res+= (0.25 * prs->cfuncderiv(nbr,x,y) * (tmp[0] + tmp[3]) / nbr) * (tmp[0] + tmp[3]) / 16.0;
									
				result->a[ctr] = res;
				ctr++;
			};

			//8
			res = 0;
			result->ja[ctr] = prs->eynum(i, j);
			//top
			x = prs->h/2.0 + i * prs->h;
			y = prs->h/2.0 + j * prs->h;


			res += prs->cfunc(ntl,x,y) / 16.0 + prs->cfunc(ntr,x,y) / 16.0 + prs->cfunc(ntc,x,y) / 4.0;
				
			res += (prs->cfuncderiv(ntc,x,y) * tmp[4] / ntc) * tmp[4] / 4.0;
			res += (0.25 * prs->cfuncderiv(ntr,x,y) * (tmp[4] + tmp[7]) / ntr) * (tmp[4] + tmp[7]) / 16.0;
			res += (0.25 * prs->cfuncderiv(ntl,x,y) * (tmp[4] + tmp[7]) / ntl) * (tmp[4] + tmp[7]) / 16.0;
			
			x = prs->h/2.0 + i * prs->h;
			y = prs->h/2.0 + (j - 1) * prs->h;
			if (y < 0) y+= 1.0;

			res += prs->cfunc(nbl,x,y) / 16.0 + prs->cfunc(nbr,x,y) / 16.0 + prs->cfunc(nbc,x,y) / 4.0;
				
			res += (prs->cfuncderiv(nbc,x,y) * tmp[0] / nbc) * tmp[0] / 4.0;
			res += (0.25 * prs->cfuncderiv(nbl,x,y) * (tmp[0] + tmp[3]) / nbl) * (tmp[0] + tmp[3]) / 16.0;
			res += (0.25 * prs->cfuncderiv(nbr,x,y) * (tmp[0] + tmp[3]) / nbr) * (tmp[4] + tmp[3]) / 16.0;
			
			result->a[ctr] = res;
			ctr++;

			//9
			if (j != prs->n - 1)
			{
				res = 0;
				result->ja[ctr] = prs->eynum(i, modL(j + 1, prs->n));
				x = prs->h/2.0 + i * prs->h;
				y = prs->h/2.0 + j * prs->h;
				if (y > 1.0) y -= 1.0;

				res+= (prs->cfunc(ntl,x,y) + prs->cfunc(ntr,x,y)) / 16.0;

				res+= (0.25 * prs->cfuncderiv(ntl,x,y) * (tmp[4] + tmp[7]) / ntl) * (tmp[4] + tmp[7]) / 16.0;
				res+= (0.25 * prs->cfuncderiv(ntr,x,y) * (tmp[4] + tmp[7]) / ntr) * (tmp[4] + tmp[7]) / 16.0;

				result->a[ctr] = res;
				ctr++;
			};

			//7, loop			
			if (j == 0)
			{
				x = prs->h/2.0 + i * prs->h;
				y = prs->h/2.0 + (j - 1) * prs->h;
				if (y < 0) y += 1.0;

				res = 0;
				result->ja[ctr] = prs->eynum(i, modL(j - 1, prs->n));

				res+= (prs->cfunc(nbl,x,y) + prs->cfunc(nbr,x,y)) / 16.0;

				res+= (0.25 * prs->cfuncderiv(nbl,x,y) * (tmp[0] + tmp[3]) / nbl) * (tmp[0] + tmp[3]) / 16.0;
				res+= (0.25 * prs->cfuncderiv(nbr,x,y) * (tmp[0] + tmp[3]) / nbr) * (tmp[0] + tmp[3]) / 16.0;
									
				result->a[ctr] = res;
				ctr++;
			};

			//gradientY
			//1
			result->a[ctr] = +1.0/prs->h;
			if (j != 0)
			{
				result->a[ctr] = +1.0 * prs->scaling/prs->h;
				result->ja[ctr] = prs->pnum(i, j - 1);
				ctr++;
			};

			//2
			result->ja[ctr] = prs->pnum(i, j);
			result->a[ctr] = -1.0 * prs->scaling / prs->h;
			ctr++;

			//1.5
			if (j == 0)
			{
				result->a[ctr] = +1.0 * prs->scaling /prs->h;
				result->ja[ctr] = prs->pnum(i, prs->n - 1);
				ctr++;
			};

			result->ia[eqnum + 1] = ctr;
			eqnum++;
		};
	};
	

	//simple equation
	for(int i = 0; i < 2 * prs->n; i++)
	{
		for(int j = 0; j < prs->n; j++)
		{
			//1
			result->ja[ctr] = prs->exnum(i, j);
			result->a[ctr] = - prs->scaling / prs->h;
			ctr++;

			//2 
			result->ja[ctr] = prs->exnum(i + 1, j);
			result->a[ctr] = prs->scaling / prs->h;
			ctr++;


			//4 b
			if (j == prs->n - 1)
			{
				result->ja[ctr] = prs->eynum(i, modL(j + 1, prs->n));
				result->a[ctr] = prs->scaling / prs->h;
				ctr++;
			};			

			//3 (a,b)
			result->ja[ctr] = prs->eynum(i, j);
			result->a[ctr] = - prs->scaling / prs->h;
			ctr++;

			//4 a
			if (j != prs->n - 1)
			{
				result->ja[ctr] = prs->eynum(i, modL(j + 1, prs->n));
				result->a[ctr] = prs->scaling / prs->h;
				ctr++;
			};
			
			result->ia[eqnum + 1] = ctr;
			eqnum++;
		};
	};


	if(ctr != nnz)
	{
		printf("(%d) ", ctr);
		printf("Fullscheme: Derivative matrix #NNZ mismatch!\n");
	};

	for(int i = 0; i < prs->varnum + 1; i++)
	{
		result->ia[i] += 1;
	};

	for(int i = 0; i < nnz; i++)
	{
		result->ja[i] += 1;
	};

	return result;
};

matrixCSR * partderivCSR(double * pt, void * params, int rank, int * ia1)
{
	solverparams * prs;
	prs = (solverparams *) params;

	matrixCSR * result;

	result    = new matrixCSR();

	int vn = prs->varnum;
	int nnz = 9 * rank;

	result->n = vn;
	result->m = rank;
	result->a = new double[nnz];
	result->ia = new int[rank + 1];
	result->ja = new int[nnz];

	int ctr = 0;
	int eqnum = 0;
	result->ia[0] = ctr;

	int i, j;
	int num;

	double x,y;

	for(int k = 0; k < rank; k++)
	{
		num = ia1[k];
		if (num < (2 * prs->n + 1) * prs->n)
		{
			i = num / prs->n;
			j = num - i * prs->n;

			//nl equation, x
			double tmp[8];
			double nlt, nlb, nlc, nrt, nrb, nrc;
			double res;
			if (i != 0)
			{
				tmp[0] = pt[prs->exnum(i, j)];
				tmp[3] = pt[prs->exnum(i - 1, j)];
		
				tmp[1] = pt[prs->eynum(i - 1, j)];
				tmp[2] = pt[prs->eynum(i - 1, modL(j + 1, prs->n))];
			};
			
			if (i != 2 * prs->n)
			{
				tmp[4] = pt[prs->exnum(i, j)];
				tmp[7] = pt[prs->exnum(i + 1, j)];
			
				tmp[5] = pt[prs->eynum(i, j)];
				tmp[6] = pt[prs->eynum(i, modL(j+1, prs->n))];
			};
			
			if (i != 0)
			{
				nlc = sqrt(tmp[0] * tmp[0] + 0.25 * (tmp[1] + tmp[2]) * (tmp[1] + tmp[2]));
				nlb = sqrt(tmp[1] * tmp[1] + 0.25 * (tmp[0] + tmp[3]) * (tmp[0] + tmp[3]));
				nlt = sqrt(tmp[2] * tmp[2] + 0.25 * (tmp[0] + tmp[3]) * (tmp[0] + tmp[3]));
			};

			if (i != 2 * prs->n)
			{
				nrc = sqrt(tmp[4] * tmp[4] + 0.25 * (tmp[5] + tmp[6]) * (tmp[5] + tmp[6]));
				nrb = sqrt(tmp[5] * tmp[5] + 0.25 * (tmp[4] + tmp[7]) * (tmp[4] + tmp[7]));
				nrt = sqrt(tmp[6] * tmp[6] + 0.25 * (tmp[4] + tmp[7]) * (tmp[4] + tmp[7]));
			};
				
			//3
			x = prs->h/2.0 + (i-1) * prs->h;
			y = prs->h/2.0 + j * prs->h;

			if (i != 0)
			{
				res = 0;
				result->ja[ctr] = prs->exnum(i-1, j);
					
				res+= (prs->cfunc(nlt,x,y) + prs->cfunc(nlb,x,y)) / 16.0;
					
				res+= (0.25 * prs->cfuncderiv(nlt,x,y) * (tmp[0] + tmp[3]) / nlt) * (tmp[0] + tmp[3]) / 16.0;
				res+= (0.25 * prs->cfuncderiv(nlb,x,y) * (tmp[0] + tmp[3]) / nlb) * (tmp[0] + tmp[3]) / 16.0;
				result->a[ctr] = res;
				ctr++;
			};

			//4
			res = 0;
			result->ja[ctr] = prs->exnum(i, j);
			//left
			x = prs->h/2.0 + (i-1) * prs->h;
			y = prs->h/2.0 + j * prs->h;

			if (i != 0)
			{
				res += prs->cfunc(nlt,x,y) / 16.0 + prs->cfunc(nlb,x,y) / 16.0 + prs->cfunc(nlc,x,y) / 4.0;
					
				res += (prs->cfuncderiv(nlc,x,y) * tmp[0] / nlc) * tmp[0] / 4.0;
				res += (0.25 * prs->cfuncderiv(nlb,x,y) * (tmp[0] + tmp[3]) / nlb) * (tmp[0] + tmp[3]) / 16.0;
				res += (0.25 * prs->cfuncderiv(nlt,x,y) * (tmp[0] + tmp[3]) / nlt) * (tmp[0] + tmp[3]) / 16.0;
			};
			
			//right
			x = prs->h/2.0 + i * prs->h;
			y = prs->h/2.0 + j * prs->h;

			if (i != 2 * prs->n)
			{
				res += prs->cfunc(nrt,x,y) / 16.0 + prs->cfunc(nrb,x,y) / 16.0 + prs->cfunc(nrc,x,y) / 4.0;
					
				res += (prs->cfuncderiv(nrc,x,y) * tmp[4] / nrc) * tmp[4] / 4.0;
				res += (0.25 * prs->cfuncderiv(nrb,x,y) * (tmp[4] + tmp[7]) / nrb) * (tmp[4] + tmp[7]) / 16.0;
				res += (0.25 * prs->cfuncderiv(nrt,x,y) * (tmp[4] + tmp[7]) / nrt) * (tmp[4] + tmp[7]) / 16.0;
			};
			result->a[ctr] = res;
			ctr++;
	
			//5
			x = prs->h/2.0 + i * prs->h;
			y = prs->h/2.0 + j * prs->h;
			if (i != 2 * prs->n)
			{
				res = 0;
				result->ja[ctr] = prs->exnum(i + 1, j);

				res+= (prs->cfunc(nrt,x,y) + prs->cfunc(nrb,x,y)) / 16.0;
				
				res+= (0.25 * prs->cfuncderiv(nrt,x,y) * (tmp[4] + tmp[7]) / nrt) * (tmp[4] + tmp[7]) / 16.0;
				res+= (0.25 * prs->cfuncderiv(nrb,x,y) * (tmp[4] + tmp[7]) / nrb) * (tmp[4] + tmp[7]) / 16.0;
				result->a[ctr] = res;
				ctr++;
			};

			//6
			x = prs->h/2.0 + (i - 1) * prs->h;
			y = prs->h/2.0 + j * prs->h;
			if (i != 0)
			{
				if (j != prs->n - 1)
				{
					res = 0;
					result->ja[ctr] = prs->eynum(i - 1, j);

					res += (prs->cfuncderiv(nlb,x,y) * tmp[1] / nlb) * (tmp[0] + tmp[3]) / 16.0;
					res += (prs->cfuncderiv(nlc,x,y) * 0.25 * (tmp[1] + tmp[2]) / nlc) * tmp[0] / 4.0;

					result->a[ctr] = res;
					ctr++;
				};
			};

			//7
			if (i != 0)
			{
				res = 0;
				result->ja[ctr] = prs->eynum(i - 1, modL(j + 1, prs->n));

				res += (prs->cfuncderiv(nlt,x,y) * tmp[2] / nlt) * (tmp[0] + tmp[3]) / 16.0;
				res += (prs->cfuncderiv(nlc,x,y) * 0.25 * (tmp[1] + tmp[2]) / nlc) * tmp[0] / 4.0;
			
				result->a[ctr] = res;
				ctr++;
			};

			//6, loop
			if (i != 0)
			{
				if (j == prs->n - 1)
				{
					res = 0;
					result->ja[ctr] = prs->eynum(i - 1, j);

					res += (prs->cfuncderiv(nlb,x,y) * tmp[1] / nlb) * (tmp[0] + tmp[3]) / 16.0;
					res += (prs->cfuncderiv(nlc,x,y) * 0.25 * (tmp[1] + tmp[2]) / nlc) * tmp[0] / 4.0;

					result->a[ctr] = res;
					ctr++;
				};
			};

			//8
			x = prs->h/2.0 + i * prs->h;
			y = prs->h/2.0 + j * prs->h;

			if (i != 2 * prs->n)
			{
				if (j != prs->n - 1)
				{
					res = 0;
					result->ja[ctr] = prs->eynum(i, j);

					res += (prs->cfuncderiv(nrb,x,y) * tmp[5] / nrb) * (tmp[4] + tmp[7]) / 16.0;
					res += (prs->cfuncderiv(nrc,x,y) * 0.25 * (tmp[5] + tmp[6]) / nrc) * tmp[4] / 4.0;

					result->a[ctr] = res;
					ctr++;
				};
			};

			//9
			if (i != 2 * prs->n)
			{
				res = 0;
				result->ja[ctr] = prs->eynum(i, modL(j + 1, prs->n));

				res += (prs->cfuncderiv(nrt,x,y) * tmp[6] / nrt) * (tmp[4] + tmp[7]) / 16.0;
				res += (prs->cfuncderiv(nrc,x,y) * 0.25 * (tmp[5] + tmp[6]) / nrc) * tmp[4] / 4.0;
				
				result->a[ctr] = res;
				ctr++;
			};

			//8, loop
			if (i != 2 * prs->n)
			{
				if (j == prs->n - 1)
				{
					res = 0;
					result->ja[ctr] = prs->eynum(i, j);

					res += (prs->cfuncderiv(nrb,x,y) * tmp[5] / nrb) * (tmp[4] + tmp[7]) / 16.0;
					res += (prs->cfuncderiv(nrc,x,y) * 0.25 * (tmp[5] + tmp[6]) / nrc) * tmp[4] / 4.0;

					result->a[ctr] = res;
					ctr++;
				};
			};

		}
		else
		{
			num -= (2 * prs->n + 1) * prs->n;
			if (num < (prs->n) * 2 * prs->n)
			{
				i = num / (prs->n);
				j = num - i * (prs->n);

				//nl equation, y

				double tmp[8];
				double ntc, ntl, ntr, nbc, nbl, nbr;
				double res;

				tmp[0] = pt[prs->eynum(i, j)];
				tmp[3] = pt[prs->eynum(i, modL(j - 1, prs->n))];

				tmp[1] = pt[prs->exnum(i, modL(j - 1, prs->n))];
				tmp[2] = pt[prs->exnum(i + 1, modL(j - 1, prs->n))];

				tmp[4] = pt[prs->eynum(i, j)];
				tmp[7] = pt[prs->eynum(i, modL(j + 1, prs->n))];

				tmp[5] = pt[prs->exnum(i, j)];
				tmp[6] = pt[prs->exnum(i + 1, j)];

				nbc = sqrt(tmp[0] * tmp[0] + 0.25 * (tmp[1] + tmp[2]) * (tmp[1] + tmp[2]));
				nbl = sqrt(tmp[1] * tmp[1] + 0.25 * (tmp[0] + tmp[3]) * (tmp[0] + tmp[3]));
				nbr = sqrt(tmp[2] * tmp[2] + 0.25 * (tmp[0] + tmp[3]) * (tmp[0] + tmp[3]));

				ntc = sqrt(tmp[4] * tmp[4] + 0.25 * (tmp[5] + tmp[6]) * (tmp[5] + tmp[6]));
				ntl = sqrt(tmp[5] * tmp[5] + 0.25 * (tmp[4] + tmp[7]) * (tmp[4] + tmp[7]));
				ntr = sqrt(tmp[6] * tmp[6] + 0.25 * (tmp[4] + tmp[7]) * (tmp[4] + tmp[7]));

				//3
				if (j != 0)
				{
					res = 0;

					x = prs->h/2.0 + i * prs->h;
					y = prs->h/2.0 + (j - 1) * prs->h;

					if (y < 0) y+= 1.0;

					result->ja[ctr] = prs->exnum(i, modL(j - 1, prs->n));

					res += (prs->cfuncderiv(nbl,x,y) * tmp[1] / nbl) * (tmp[0] + tmp[3]) / 16.0;
					res += (prs->cfuncderiv(nbc,x,y) * 0.25 * (tmp[1] + tmp[2]) / nbc) * tmp[0] / 4.0;

					result->a[ctr] = res;
					ctr++;
				};

				//4
				x = prs->h/2.0 + i * prs->h;
				y = prs->h/2.0 + j * prs->h;

				res = 0;
				result->ja[ctr] = prs->exnum(i, j);

				res += (prs->cfuncderiv(ntl,x,y) * tmp[5] / ntl) * (tmp[4] + tmp[7]) / 16.0;
				res += (prs->cfuncderiv(ntc,x,y) * 0.25 * (tmp[5] + tmp[6]) / ntc) * tmp[4] / 4.0;

				result->a[ctr] = res;
				ctr++;

				//3.5
				if (j == 0)
				{
					x = prs->h/2.0 + i * prs->h;
					y = prs->h/2.0 + (j - 1) * prs->h;

					if (y < 0) y+= 1.0;
					res = 0;

					result->ja[ctr] = prs->exnum(i, modL(j-1, prs->n));

					res += (prs->cfuncderiv(nbl,x,y) * tmp[1] / nbl) * (tmp[0] + tmp[3]) / 16.0;
					res += (prs->cfuncderiv(nbc,x,y) * 0.25 * (tmp[1] + tmp[2]) / nbc) * tmp[0] / 4.0;

					result->a[ctr] = res;
					ctr++;
				};


				//5
				if (j != 0)
				{
					x = prs->h/2.0 + i * prs->h;
					y = prs->h/2.0 + (j-1) * prs->h;
					if (y < 0) y+=1.0;

					res =  0;
					result->ja[ctr] = prs->exnum(i+1, modL(j - 1, prs->n));

					res += (prs->cfuncderiv(nbr,x,y) * tmp[2] / nbr) * (tmp[0] + tmp[3]) / 16.0;
					res += (prs->cfuncderiv(nbc,x,y) * 0.25 * (tmp[1] + tmp[2]) / nbc) * tmp[0] / 4.0;

					result->a[ctr] = res;
					ctr++;
				};

				//6
				res = 0;
				result->ja[ctr] = prs->exnum(i + 1, j);
				x = prs->h/2.0 + i * prs->h;
				y = prs->h/2.0 + j * prs->h;

				res += (prs->cfuncderiv(ntr,x,y) * tmp[6] / ntr) * (tmp[4] + tmp[7]) / 16.0;
				res += (prs->cfuncderiv(ntc,x,y) * 0.25 * (tmp[5] + tmp[6]) / ntc) * tmp[4] / 4.0;


				result->a[ctr] = res;
				ctr++;

				//5.5
				if (j == 0)
				{
					res =  0;

					x = prs->h/2.0 + i * prs->h;
					y = prs->h/2.0 + (j-1) * prs->h;

					if (y < 0) y+=1.0;
					result->ja[ctr] = prs->exnum(i+1, modL(j - 1, prs->n));

					res += (prs->cfuncderiv(nbr,x,y) * tmp[2] / nbr) * (tmp[0] + tmp[3]) / 16.0;
					res += (prs->cfuncderiv(nbc,x,y) * 0.25 * (tmp[1] + tmp[2]) / nbc) * tmp[0] / 4.0;

					result->a[ctr] = res;
					ctr++;
				};

				//9, loop
				if (j == prs->n - 1)
				{
					res = 0;
					result->ja[ctr] = prs->eynum(i, modL(j + 1, prs->n));
					x = prs->h/2.0 + i * prs->h;
					y = prs->h/2.0 + j * prs->h;
					if (y > 1.0) y -= 1.0;

					res+= (prs->cfunc(ntl,x,y) + prs->cfunc(ntr,x,y)) / 16.0;

					res+= (0.25 * prs->cfuncderiv(ntl,x,y) * (tmp[4] + tmp[7]) / ntl) * (tmp[4] + tmp[7]) / 16.0;
					res+= (0.25 * prs->cfuncderiv(ntr,x,y) * (tmp[4] + tmp[7]) / ntr) * (tmp[4] + tmp[7]) / 16.0;

					result->a[ctr] = res;
					ctr++;
				};

				//7			
				if (j != 0)
				{
					x = prs->h/2.0 + i * prs->h;
					y = prs->h/2.0 + (j - 1) * prs->h;
					if (y < 0) y += 1.0;

					res = 0;
					result->ja[ctr] = prs->eynum(i, modL(j - 1, prs->n));

					res+= (prs->cfunc(nbl,x,y) + prs->cfunc(nbr,x,y)) / 16.0;

					res+= (0.25 * prs->cfuncderiv(nbl,x,y) * (tmp[0] + tmp[3]) / nbl) * (tmp[0] + tmp[3]) / 16.0;
					res+= (0.25 * prs->cfuncderiv(nbr,x,y) * (tmp[0] + tmp[3]) / nbr) * (tmp[0] + tmp[3]) / 16.0;

					result->a[ctr] = res;
					ctr++;
				};

				//8
				res = 0;
				result->ja[ctr] = prs->eynum(i, j);
				//top
				x = prs->h/2.0 + i * prs->h;
				y = prs->h/2.0 + j * prs->h;


				res += prs->cfunc(ntl,x,y) / 16.0 + prs->cfunc(ntr,x,y) / 16.0 + prs->cfunc(ntc,x,y) / 4.0;

				res += (prs->cfuncderiv(ntc,x,y) * tmp[4] / ntc) * tmp[4] / 4.0;
				res += (0.25 * prs->cfuncderiv(ntr,x,y) * (tmp[4] + tmp[7]) / ntr) * (tmp[4] + tmp[7]) / 16.0;
				res += (0.25 * prs->cfuncderiv(ntl,x,y) * (tmp[4] + tmp[7]) / ntl) * (tmp[4] + tmp[7]) / 16.0;

				x = prs->h/2.0 + i * prs->h;
				y = prs->h/2.0 + (j - 1) * prs->h;
				if (y < 0) y+= 1.0;

				res += prs->cfunc(nbl,x,y) / 16.0 + prs->cfunc(nbr,x,y) / 16.0 + prs->cfunc(nbc,x,y) / 4.0;

				res += (prs->cfuncderiv(nbc,x,y) * tmp[0] / nbc) * tmp[0] / 4.0;
				res += (0.25 * prs->cfuncderiv(nbl,x,y) * (tmp[0] + tmp[3]) / nbl) * (tmp[0] + tmp[3]) / 16.0;
				res += (0.25 * prs->cfuncderiv(nbr,x,y) * (tmp[0] + tmp[3]) / nbr) * (tmp[4] + tmp[3]) / 16.0;

				result->a[ctr] = res;
				ctr++;

				//9
				if (j != prs->n - 1)
				{
					res = 0;
					result->ja[ctr] = prs->eynum(i, modL(j + 1, prs->n));
					x = prs->h/2.0 + i * prs->h;
					y = prs->h/2.0 + j * prs->h;
					if (y > 1.0) y -= 1.0;

					res+= (prs->cfunc(ntl,x,y) + prs->cfunc(ntr,x,y)) / 16.0;

					res+= (0.25 * prs->cfuncderiv(ntl,x,y) * (tmp[4] + tmp[7]) / ntl) * (tmp[4] + tmp[7]) / 16.0;
					res+= (0.25 * prs->cfuncderiv(ntr,x,y) * (tmp[4] + tmp[7]) / ntr) * (tmp[4] + tmp[7]) / 16.0;

					result->a[ctr] = res;
					ctr++;
				};

				//7, loop			
				if (j == 0)
				{
					x = prs->h/2.0 + i * prs->h;
					y = prs->h/2.0 + (j - 1) * prs->h;
					if (y < 0) y += 1.0;

					res = 0;
					result->ja[ctr] = prs->eynum(i, modL(j - 1, prs->n));

					res+= (prs->cfunc(nbl,x,y) + prs->cfunc(nbr,x,y)) / 16.0;

					res+= (0.25 * prs->cfuncderiv(nbl,x,y) * (tmp[0] + tmp[3]) / nbl) * (tmp[0] + tmp[3]) / 16.0;
					res+= (0.25 * prs->cfuncderiv(nbr,x,y) * (tmp[0] + tmp[3]) / nbr) * (tmp[0] + tmp[3]) / 16.0;

					result->a[ctr] = res;
					ctr++;
				};
			}
			else
			{
				num -= 2 * prs->n * (prs->n);

				i = num / prs->n;
				j = num - i * prs->n;

				//simple equation
				// nothing here =)
			};
		};

		result->ia[eqnum + 1] = ctr;
		eqnum++;
	};

	int nnzreal = result->ia[rank];
	for(int i = 0; i < rank + 1; i++)
	{
		result->ia[i] += 1;
	};

	for(int i = 0; i < nnzreal; i++)
	{
		result->ja[i] += 1;
	};


	return result;
};

/*
int main(int argc, char ** argv)
{
};*/


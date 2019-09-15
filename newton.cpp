#include <mkl.h>
#include <mkl_rci.h>
#include <mkl_blas.h>
#include <mkl_spblas.h>
#include <mkl_service.h>
#include <cmath>
#include <cstdio>
#include "matrixCSR.h"
#include "newton.h"
#include "constants.h"

void newton(int n, double * point, eval function, derivCSR derivmaker, void * params, double tolerance, int maxit, double startsigma, double lstolerance, int method, int covmethod)
{
	int ione = 1;
	int izero = 0;
	double dzero = 0.0;
	double dminusone = -1.0;

	int dummy;

	double * ev = new double[n];	
	double * snu = new double[n];
	double * snuinv = new double[n];
	double * mtev      = new double[n];
	double * evtest    = new double[n];
	double * pointtest = new double[n];
        double * pointincsave = new double[n];
        double * evincsave    = new double[n];

	double NU, SIGMA;

        double residual, startresidual, newresidual, residualincsave;
        double angle, stepnorm, nuchange, scalinv, newsigma;

	bool acceptable, increase, truenewton;

        int iteration = 0;
        int ierr;

	matrixCSR * deriv = NULL;

	function(point, ev, params);
	residual = dnrm2(&n, ev, &ione);
	startresidual = residual;
	deriv = derivmaker(point, params);
	SIGMA = startsigma * startresidual;

//nustarters = residual;
	//main newton iteration cycle.
	while((residual > tolerance * startresidual) && (iteration < maxit))
	{
		//deriv->print();
		//scanf("%d", &dummy);
	
		acceptable = false;
		increase   = false;

                while(!acceptable)
                {
			truenewton = true;
                        /* SEARHING FOR (S(0)) */

                        //solving the system

                        if (method == METHOD_PRECGMRES)
                        {
                                deriv->solveprecgmres(snu, ev, lstolerance);
                        }
                        else
                        {
                                if (method = METHOD_SVD)
                                {
                                        deriv->solvesvd(snu, ev);
                                }
                                else
                                {
                                        printf("Wrong method specified for newton.\n");
                                };
                        };

                        stepnorm = dnrm2(&n, snu, &ione);
			deriv->leftmatvec(ev, mtev);

                        if (stepnorm > SIGMA)
                        {
				truenewton = false;
                                //deriv->leftmatvec(ev, mtev);
                                NU = 0;

                                /*SEARCHING FOR A NU THAT CORRESPONDS TO THE CURRENT SIGMA */
                                do
                                {
                                        scalinv = 1.0 / stepnorm;

                                        dscal(&n, &scalinv, snu, &ione);

                                        /* SOLVING NOISY HESSIAN SYSTEM WITH SNU */					
                                       
					if (covmethod == METHOD_COVSVD)
					{
						deriv->solvecovsvd(snuinv, snu, NU);
					}
					else
					{
						if (covmethod == METHOD_PRECCOVCG)
						{
							deriv->solvepreccovcg(snuinv, snu, NU, lstolerance);
						}
						else
						{
							printf("Wrong method specified for newton.\n");
						};
					};

					/*UPDATING NU */
                                        nuchange = 1.0 / ddot(&n, snu, &ione, snuinv, &ione);
                                        if (printlevel >= 2) printf(" Rayleigh: %e", ddot(&n, snu, &ione, snuinv, &ione));
                                        nuchange *= (stepnorm - SIGMA) / SIGMA;
                                        if (printlevel >= 2) printf("#NU: %e ->", NU);
                                        NU += nuchange;
                                        if (printlevel >= 2) printf("%e\n", NU);

                                        /* SEARHING FOR (S(NU)) */

					if (covmethod == METHOD_COVSVD)
					{
						deriv->solvecovsvd(snu, mtev, NU);
					}
					else
					{
						if (covmethod == METHOD_PRECCOVCG)
						{
							deriv->solvepreccovcg(snu, mtev, NU, lstolerance);
						}
						else
						{
							printf("Wrong method specified for newton.\n");
						};
					};

                                        stepnorm = dnrm2(&n, snu, &ione);
                                }
                                while ((stepnorm < 0.75 * SIGMA) || (stepnorm > 1.5 * SIGMA));
                        };

			//a point in trust region found, validating trust region
		
			//angle = ddot(&n, ev, &ione, snu, &ione);
		
			dcopy(&n, point, &ione, pointtest, &ione);
			daxpy(&n, &dminusone, snu, &ione, pointtest, &ione);
			
			function(pointtest, evtest, params);
			newresidual = dnrm2(&n, evtest, &ione);

			if (newresidual < residual)
			{
				//good point, can continue or even expand the region

                                //checking if current sigma can be increased

				double predicted, actual;

				//using snuinv as temporary variable here
				deriv->matvec(snu, snuinv);
				daxpy(&n, &dminusone, ev, &ione, snuinv, &ione);

				actual    = newresidual * newresidual - residual * residual;

				predicted = ddot(&n, snuinv, &ione, snuinv, &ione) - residual * residual;

                                if (printlevel >= 2) printf("(%e %e)\n", actual, predicted);
				if ((std::abs(predicted - actual) < 0.1 * std::abs(actual)) && (truenewton == false))
				{
					//trying to increase sigma
					increase = true;
                                        dcopy(&n, pointtest, &ione, pointincsave, &ione);
                                        dcopy(&n, evtest, &ione, evincsave, &ione);

                                        residualincsave = newresidual;
                                        
                                        SIGMA *= 2.0;
                                        if (printlevel >= 2) printf("#Sigma increase attempt: %e\n", SIGMA);
				}
				else
				{
					//good point, but no trust region increase
					acceptable = true;
                                        dcopy(&n, pointtest, &ione, point, &ione);
                                        dcopy(&n, evtest, &ione, ev, &ione);

                                        delete deriv;
                                        deriv = derivmaker(point, params);

					 // deriv->print();
					 // scanf("%d", &dummy);

					if (truenewton)
					{
						SIGMA = stepnorm;
						if (printlevel >= 2) printf("#Sigma decreased after true newton step: %e\n", SIGMA);
					};

                                        residual = newresidual;
				};
			}
			else
			{
				if (increase)
				{
					//attempt to increase region failed
					increase = false;
					SIGMA /= 2.0;

                                        if (printlevel >= 2) printf("#Sigma increase attempt failed: %e\n", SIGMA);

                                        residual = residualincsave;
                                        dcopy(&n, pointincsave, &ione, point, &ione);
                                        dcopy(&n, evincsave, &ione, ev, &ione);

                                        delete deriv;
                                        deriv = derivmaker(point, params);

					 // deriv->print();
					 // scanf("%d", &dummy);
                                        acceptable = true;
				}
				else
				{
					//bad point, shrinking trust region
					double lambda;
					angle = - ddot(&n, mtev, &ione, snu, &ione);

					if (newresidual * newresidual - residual * residual - 2.0 * angle != 0.0)
					{
						lambda = - angle / (newresidual * newresidual - residual * residual - 2.0 * angle);
						newsigma = lambda * stepnorm; //dnrm2(&n, snu, &ione);
					}
					else
					{
						newsigma = SIGMA;	
					};

					if (newsigma < 0.1 * SIGMA) newsigma = 0.1 * SIGMA;
					if (newsigma > 0.5 * SIGMA) newsigma = 0.5 * SIGMA;

					SIGMA = newsigma;
					if (printlevel >= 2) printf("#Sigma shrinked: %e\n", SIGMA);
				
				};
			};	
			if (SIGMA < 1.0e-14 * startsigma) break;
		};

		if (SIGMA < 1.0e-14 * startsigma)
		{
			if (printlevel >= 1) printf("Newton: SIGMA converged to zero (appears to be a local minimum).\n");
			break;
		};

		if (printlevel >= 1) printf("Newton: After step %d: %e\n", iteration + 1, residual / startresidual);
		iteration++;
	};

	delete deriv;

	delete[] snu;
	delete[] snuinv;
	delete[] ev;
	delete[] evtest;
	delete[] pointtest;
        delete[] mtev;
        delete[] evincsave;
        delete[] pointincsave;
};

void newtonmark2(int n, double * point, eval function, derivCSR derivmaker, void * params, double tolerance, int maxit, double startsigma, double lstolerance, int method)
{
	int ione = 1;
	int izero = 0;
	double dzero = 0.0;
	double dminusone = -1.0;

	int dummy;

	double * ev = new double[n];	
	double * snu = new double[n];
	double * temp = new double[n];
	double * temp2 = new double[n];
	double * mtev      = new double[n];
	double * evtest    = new double[n];
	double * pointtest = new double[n];
        double * pointincsave = new double[n];
        double * evincsave    = new double[n];
	
	double SIGMA;

        double residual, startresidual, newresidual, residualincsave;
        double angle, stepnorm, nuchange, scalinv, newsigma;
	double forwardhess, backwardhess;

	bool acceptable, increase, truenewton;

        int iteration = 0;
        int ierr;

	matrixCSR * deriv = NULL;

	function(point, ev, params);
	residual = dnrm2(&n, ev, &ione);
	startresidual = residual;
	deriv = derivmaker(point, params);
	SIGMA = startsigma * startresidual;

					 // deriv->print();
					 // scanf("%d", &dummy);
	//nustarters = residual;
	//main newton iteration cycle.
	while((residual > tolerance * startresidual) && (iteration < maxit))
	{
		int dm;
		//deriv->print();
		//scanf("%d\n", &dm);
		acceptable = false;
		increase   = false;

                while(!acceptable)
                {
			truenewton = true;
                        /* SEARHING FOR (S(0)) */

                        //solving the system

                        if (method == METHOD_PRECGMRES)
                        {
                                deriv->solveprecgmres(snu, ev, lstolerance);
                        }
                        else
                        {
                                if (method = METHOD_SVD)
                                {
                                        deriv->solvesvd(snu, ev);
                                }
                                else
                                {
                                        printf("Wrong method specified for newton.\n");
                                };
                        };

                        stepnorm = dnrm2(&n, snu, &ione);
			deriv->leftmatvec(ev,mtev);
				
                        if (stepnorm > SIGMA)
                        {
				truenewton = false;
				scalinv = 1.0 / dnrm2(&n, mtev, &ione);
				dscal(&n, &scalinv, mtev, &ione);
					
				deriv->matvec(mtev, temp);
				deriv->leftmatvec(temp, temp2);
				forwardhess = 1.0 / ddot(&n, mtev, &ione, temp2, &ione);

				if (forwardhess / scalinv >= SIGMA)
				{
					//printf("???\n\n");
					dcopy(&n, mtev, &ione, snu, &ione);
					dscal(&n, &SIGMA, snu, &ione);
				}
				else
				{
					double etta;
					backwardhess = 1.0 / (scalinv * ddot(&n, mtev, &ione, snu, &ione));
					etta = 0.2 + 0.8 * backwardhess * forwardhess;
	
					//etta = 0.01 + 0.99 * backwardhess * forwardhess;
					
					if (printlevel >= 2) printf("%e|%e|%e|%e|%e\n", SIGMA, forwardhess / scalinv, forwardhess, backwardhess, etta);	
						
					if (etta * stepnorm > SIGMA)
					{		
						//printf("%e|%e|%e|%e|%e\n", SIGMA, forwardhess / scalinv, forwardhess, backwardhess, etta);	
						//temp = lambda * mtev, vec1
						scalinv = forwardhess / scalinv;
						dcopy(&n, mtev, &ione, temp, &ione);
						dscal(&n, &scalinv, temp, &ione);

						//temp2 = etta * snu - vec1
						dcopy(&n, temp, &ione, temp2, &ione);
						dscal(&n, &dminusone, temp2, &ione);
						daxpy(&n, &etta, snu, &ione, temp2, &ione);

						double a,b,c;
						a = dnrm2(&n, temp2, &ione);
						a = a * a;
						b = 2 * ddot(&n, temp, &ione, temp2, &ione);
						c = dnrm2(&n, temp, &ione);
						c = c * c - SIGMA * SIGMA;

						double stp;
						stp = (-b + sqrt(b * b - 4.0 * a * c)) / (2.0 * a);

						dcopy(&n, temp, &ione, snu, &ione);
						daxpy(&n, &stp, temp2, &ione, snu, &ione);
					}
					else
					{
						double sc  = SIGMA / stepnorm;
						if (printlevel >= 2) printf("#N#");	
						dscal(&n, &sc, snu, &ione);
					};
                            	};
	                };

			//a point in trust region found, validating trust region
		
			//angle = ddot(&n, ev, &ione, snu, &ione);
		
			dcopy(&n, point, &ione, pointtest, &ione);
			daxpy(&n, &dminusone, snu, &ione, pointtest, &ione);
			
			function(pointtest, evtest, params);
			newresidual = dnrm2(&n, evtest, &ione);

			if (newresidual < residual)
			{
				//good point, can continue or even expand the region

                                //checking if current sigma can be increased

				double predicted, actual;

				deriv->matvec(snu, temp);
				daxpy(&n, &dminusone, ev, &ione, temp, &ione);

				actual    = newresidual * newresidual - residual * residual;

				predicted = ddot(&n, temp, &ione, temp, &ione) - residual * residual;

                                if (printlevel >= 2) printf("(%e %e)\n", actual, predicted);
				if ((std::abs(predicted - actual) < 0.1 * std::abs(actual)) && (truenewton == false))
				{
					//trying to increase sigma
					increase = true;
                                        dcopy(&n, pointtest, &ione, pointincsave, &ione);
                                        dcopy(&n, evtest, &ione, evincsave, &ione);

                                        residualincsave = newresidual;
                                        
                                        SIGMA *= 2.0;
                                        if (printlevel >= 2) printf("#Sigma increase attempt: %e\n", SIGMA);
				}
				else
				{
					//good point, but no trust region increase
					acceptable = true;
                                        dcopy(&n, pointtest, &ione, point, &ione);
                                        dcopy(&n, evtest, &ione, ev, &ione);

                                        delete deriv;
                                        deriv = derivmaker(point, params);

					 //deriv->print();
					 //scanf("%d", &dummy);

					if (truenewton)
					{
						SIGMA = stepnorm;
						if (printlevel >= 2) printf("#Sigma decreased after true newton step: %e\n", SIGMA);
					};

                                        residual = newresidual;
				};
			}
			else
			{
				if (increase)
				{
					//attempt to increase region failed
					increase = false;
					SIGMA /= 2.0;

                                        if (printlevel >= 2) printf("#Sigma increase attempt failed: %e\n", SIGMA);

                                        residual = residualincsave;
                                        dcopy(&n, pointincsave, &ione, point, &ione);
                                        dcopy(&n, evincsave, &ione, ev, &ione);

                                        delete deriv;
                                        deriv = derivmaker(point, params);

					 // deriv->print();
					 // scanf("%d", &dummy);
                                        acceptable = true;
				}
				else
				{
					//bad point, shrinking trust region
					double lambda;
					angle = - ddot(&n, mtev, &ione, snu, &ione);
					if (newresidual * newresidual - residual * residual - 2.0 * angle != 0.0)
					{
						lambda = - angle / (newresidual * newresidual - residual * residual - 2.0 * angle);
						newsigma = lambda * stepnorm; //dnrm2(&n, snu, &ione);
					}
					else
					{
						newsigma = SIGMA;	
					};

					if (newsigma < 0.1 * SIGMA) newsigma = 0.1 * SIGMA;
					if (newsigma > 0.5 * SIGMA) newsigma = 0.5 * SIGMA;

					SIGMA = newsigma;
					if (printlevel >= 2) printf("#Sigma shrinked: %e\n", SIGMA);
				};
			};
			if (SIGMA < 1.0e-14 * startsigma) break;
		};

		if (SIGMA < 1.0e-14 * startsigma)
		{
			if (printlevel >= 1) printf("Newton: SIGMA converged to zero (appears to be a local minimum).\n");
			break;
		};

		if (printlevel >= 1) printf("Newton: After step %d: %e\n", iteration + 1, residual / startresidual);
		iteration++;
	};

	delete deriv;

	delete[] snu;
	delete[] temp;
	delete[] temp2;
	delete[] ev;
	delete[] evtest;
	delete[] pointtest;
        delete[] mtev;
        delete[] evincsave;
        delete[] pointincsave;
};


//This code file contains various geometry functions and material property functions.

#include <cmath>
#include <cstdio>
#include "constants.h"

double simpledatars(double x, double y, double alpha, double beta, double shift)
{
	double res = 0;
	double ys;
	ys = y;
	if (x > 1)
	{
		ys += shift;
		if (ys > 1.0) ys -= 1;
		
	};
	if ((ys >= 0.25) && (ys <= 0.75))
	{
		if ((x >= 1.3) && (x <= 1.7))
		{
			res = beta;
		};
		
		if ((x >= alpha) && (x <= alpha + 0.4))
		{
			res = 1.0;
		};
	};	
	return res;
};

double sqnl(double s, double gamma)
{
	double res;
	if (gamma != 0)
	{
		if (s != 0)
		{
			res = (sqrt(4.0 * gamma * s + 1) - 1.0) / (2.0 * gamma * s);
		}
		else
		{
			res = 1.0;
		};
	}
	else
	{
		return 1.0;
	};
	return res;
};

double sqnlderiv(double s, double gamma)
{
	double res, sq;
	if (gamma != 0)
	{
		if (s != 0)
		{
			sq = sqrt(4.0 * gamma * s + 1);
			res = (sq - 2 * gamma * s - 1) / (2.0 * gamma * s * s * sq);
		}
		else
		{
			res = -gamma;
		};
	}
	else
	{
		return 0;
	};
	
	return res;
};

double expsolution(double x, double y, double alpha, double beta)
{
	double res;
	res = 1.0;
	res *= exp(1.0 / (alpha * x * (x - 2.0)) + 1.0 / alpha);
	//res*= exp(1.0 / (alpha * x * (x - 2.0) * (x + 0.1)) + 1.0 / (alpha * (8569 + 421 * sqrt(421)) / 13500));
	res *= exp(1.0 / (beta * y * (y - 1.0)) + 4.0 / beta);
	return res;
};

double sinxsolution(double x, double gammal, double gammar)
{
	double res;
	double alpha, beta;
	double gamma;
	if (x < 1.5)
	{
		gamma = gammal;
	}
	else
	{
		gamma = gammar;
	};
	alpha = sqrt(2.0) * (1.0 - sqrt(1.0 + 4.0 * gamma)) / (gamma * M_PI);
	beta = (1.0 - sqrt(1.0 + 4.0 * gamma)) / (gamma * M_PI) - 1;
	//printf("%f %f %f\n", alpha, beta, gamma);
	res = alpha * sin(M_PI * x / 2.0) + beta * sin(M_PI * x); 
	return res;
};

double sinxrs(double x, double gammal, double gammar)
{
	double res;
	double alpha, beta;
	double gamma;
	double eks, ekss;
	if (x < 1.5)
	{
		gamma = gammal;
	}
	else
	{
		gamma = gammar;
	};
	alpha = sqrt(2.0) * (1.0 - sqrt(1.0 + 4.0 * gamma)) / (gamma * M_PI);
	beta = (1.0 - sqrt(1.0 + 4.0 * gamma)) / (gamma * M_PI) - 1;

	//printf("%f %f %f\n", alpha, beta, gamma);
	eks = 0.5 * M_PI * alpha * cos(M_PI * x * 0.5) + M_PI * beta * cos(M_PI * x);
	ekss = - M_PI * M_PI * (0.25 * alpha * sin(M_PI * x / 2.0) + beta * sin(M_PI * x)); 
	
	if (eks >= 0)
	{
		res = 1.0 * ekss + 2.0 * gamma * eks * ekss;
	}
	else
	{
		res = 1.0 * ekss - 2.0 * gamma * eks * ekss;
	};
	
	return res;
};

double sinxfluxx(double x, double gammal, double gammar)
{
	double res;
	
	double alpha, beta;
	double gamma;
	if (x < 1.5)
	{
		gamma = gammal;
	}
	else
	{
		gamma = gammar;
	};
	alpha = sqrt(2.0) * (1.0 - sqrt(1.0 + 4.0 * gamma)) / (gamma * M_PI);
	beta = (1.0 - sqrt(1.0 + 4.0 * gamma)) / (gamma * M_PI) - 1;
	
	res = M_PI * 0.5 * alpha * cos(M_PI * x / 2.0) + M_PI * beta * cos(M_PI * x);

	return (1.0 + gamma * std::abs(res)) * res;
};

double expX(double x, double alpha)
{
	double res;
	if (x != 0.0 && x != 2.0)
	{
		res = - 2.0 * (x - 1) / (alpha * (x - 2) * (x - 2) * x * x);
		//res = -3.0 * (x * x - 19.0 * x / 15.0 - 1.0 / 15.0) / (alpha * x * x * (x - 2) * (x - 2) * (x + 0.1) * (x + 0.1));
	}
	else
	{
		res = 0.0;
	};
	return res;
};

double expY(double y, double beta)
{
	double res;
	if (y != 0.0 && y != 1.0)
	{
		res = -2.0 * (y - 0.5) / (beta * (y - 1) * (y - 1) * y * y);
	}
	else
	{
		res = 0.0;
	};
	return res;
};

double expXs(double x,double alpha)
{
	double res;
	if (x != 0.0 && x != 2.0)
	{
		res = 2.0 * (3.0 * x * x - 6.0 * x + 4) / (alpha * (x - 2) * (x - 2) * (x - 2) * x * x * x);
		//res = 12.0 * (x * x * x * x - 38.0 * x * x * x / 15.0 + 1.705 * x * x + 0.19 * x + 1.0 / 150.0) / (alpha * (x-2) * (x-2) * (x-2) * x * x * x * (x+0.1) * (x+0.1) * (x+0.1));
	}
	else
	{
		res = 0.0;
	};
	return res;
};

double expYs(double y, double beta)
{
	double res;
	if (y != 0.0 && y != 1.0)
	{
		res = (6.0 * y * y - 6.0 * y + 2.0) / (beta * (y - 1) * (y - 1) * (y - 1) * y * y * y);
	}
	else
	{
		res = 0.0;
	};
	return res;
};

double expfluxx(double x, double y, double alpha, double beta, double gammal, double gammag, double gammar, double shift, double gap)
{
	double yp;
	double ex, ey, exs, eys, sq, ef;
	double res;
        double gamma;

        if (x < 1.0)
        {
                gamma = gammal;
        }
        else
        {
                if (x < gap)
                {
                        gamma = gammag;
                }
                else
                {
                        gamma = gammar;
                };
        };

	yp = y;
	/*
	if (x > 1.0)
	{
		if (y - shift < 0.0) 
		{
			yp = y - shift + 1.0;
		}
		else
		{
			yp = y - shift;
		};
	}
	else
	{
		yp = y;
	};*/
	ef = expsolution(x, yp, alpha, beta);
	ex = expX(x, alpha);
	ey = expY(yp, beta);
	exs = expXs(x, alpha);
	eys = expYs(yp, beta);
	sq = sqrt(ex * ex + ey * ey);
	
	res = (1.0 + gamma * ef * sq) * ef * ex;

	return res;
};

double expfluxy(double x, double y, double alpha, double beta, double gammal, double gammag, double gammar, double shift, double gap)
{
	double yp;
	double ex, ey, exs, eys, sq, ef;
	double res;
        double gamma;

        if (x < 1.0)
        {
                gamma = gammal;
        }
        else
        {
                if (x < gap)
                {
                        gamma = gammag;
                }
                else
                {
                        gamma = gammar;
                };
        };
	yp = y;
	/*
	if (x > 1.0)
	{
		if (y - shift < 0.0) 
		{
			yp = y - shift + 1.0;
		}
		else
		{
			yp = y - shift;
		};
	}
	else
	{
		yp = y;
	};*/
	ef = expsolution(x, yp, alpha, beta);
	ex = expX(x, alpha);
	ey = expY(yp, beta);
	exs = expXs(x, alpha);
	eys = expYs(yp, beta);
	sq = sqrt(ex * ex + ey * ey);
	
	res = (1.0 + gamma * ef * sq) * ef * ey;

	return res;
};

double exprs(double x, double y, double alpha, double beta, double gammal, double gammag, double gammar, double shift, double gap)
{
	double yp;
	double ex, ey, exs, eys, sq, ef;
	double res;
        double gamma;

        if (x < 1.0)
        {
                gamma = gammal;
        }
        else
        {
                if (x < gap)
                {
                        gamma = gammag;
                }
                else
                {
                        gamma = gammar;
                };
        };

	if (x > 1.0)
	{
		if (y - shift < 0.0) 
		{
			yp = y - shift + 1.0;
		}
		else
		{
			yp = y - shift;
		};
	}
	else
	{
		yp = y;
	};

	ef = expsolution(x, yp, alpha, beta);
	ex = expX(x, alpha);
	ey = expY(yp, beta);
	exs = expXs(x, alpha);
	eys = expYs(yp, beta);
	sq = sqrt(ex * ex + ey * ey);
	
	res = 0.0;
	res += (1.0 + gamma * ef * sq) * (ef * ex * ex + ef * exs);
	res += ef * ex * (gamma * ef * ex * sq + gamma * ef * ex * exs / sq);

	res += (1.0 + gamma * ef * sq) * (ef * ey * ey + ef *  eys);
	res += ef * ey * (gamma * ef * ey * sq + gamma * ef * ey *eys / sq);
	
	return res;
};

double VoidNL()
{
	return 4 * M_PI * 1.0e-7;;
};

double ArmNL(double s)
{
	/*
        double res;
        if (s != 0)
        {
                res = (0.5 * log(1 + 0.082 * s) + 4 * M_PI * 10e-7 * s) / s;
		return res;
        }
        else
        {
                return 0.410126;
        };
	*/

	double huge = 62500000000.0;
	
	if (s != 0)
	{
		return VoidNL() + 1.6 * s * s * s / (huge + s * s * s * s);
	}
	else
	{
		return VoidNL();
	};
	//return VoidNL() + 300 * s * s / (s * s * s * s + 1);
	
	//return 10 * exp(- 5 * s) + 1;
};

double ArmDeriv(double s)
{
	double tmp;
	double huge = 62500000000.0;
	/*
        double res;

        if (s != 0)
        {
		
                res = 0.5 * 0.082 / (1 + 0.082 * s) + 4 * M_PI * 10e-7;
                res -= ArmNL(s);
                res /= s;
		
		return res;
        }
        else
        {
                printf("Error: ArmDeriv referenced at zero.\n");
                return - 500000;
        };
	*/
	
	if (s != 0)
	{
		tmp = 1.6 * s * s;
		tmp /= (s * s * s * s + huge);
		tmp *= (huge * 3.0 - s * s * s * s);
		tmp /= (s * s * s * s + huge);
		return tmp;
	}
	else
	{
		return 0;
	};
	//return -50 * exp(-5 * s);
	//return 300 * 2.0 * s * (1.0 - s * s * s * s) / ((1.0 + s * s * s * s) * (1.0 + s * s * s * s));
	
};

double AirNL(double s)
{
	//return 4 * M_PI * 10e-7;
	return VoidNL();
};

double AirDeriv(double s)
{
	return 0;
};

double WoodNL(double s)
{
	return 2.0 * VoidNL();
};

double WoodDeriv(double s)
{
	return 0;
};

double IronNL(double s)
{
	/*
	double res;

	double vzero = 0.041037123051898846;

	if (s != 0)
	{
		res = 1.4 / (1.0 + exp(0.01 * (350 - s))) + log(1 + 0.00007 * s) / 12.0 + 4 * M_PI * 10e-7 * s - vzero; 
		res /= s;
		return res;
	}
	else
	{
		res = 0.000405432;
	};
	return res;
	*/

	double huge = 62500000000.0;
	
	if (s != 0)
	{
		return VoidNL() + s * s * s / (huge + s * s * s * s);
	}
	else
	{
		return VoidNL();
	};

	//return VoidNL() + 500 * s * s / (s * s * s * s + 1);
	//return 50 * exp(- 5 * s) + 1;
};

double IronDeriv(double s)
{
	double tmp;
	double huge = 62500000000.0;
	/*
	double res;
	if (s != 0)
	{
		res  = 0;
		res += 0.014 * exp(0.01 * (350 - s)) / ((1 + exp(0.01 * (350 - s))) * (1 + exp(0.01 * (350 - s))));
		res += 0.00007 / (12.0 * (1 + 0.00007 * s));
		res += 4.0 * M_PI * 10e-7;

		res -= IronNL(s);
		res /= s;
                return res;
	}
	else
	{
		printf("Error: referenced Iron Deriv at zero point.\n");
		return -5000000;
	};
	*/
	if (s != 0)
	{
		tmp = s * s;
		tmp /= (s * s * s * s + huge);
		tmp *= (huge * 3.0 - s * s * s * s);
		tmp /= (s * s * s * s + huge);
		return tmp;
	}
	else
	{
		return 0;
	};

	//return 500 * 2.0 * s * (1.0 - s * s * s * s) / ((1.0 + s * s * s * s) * (1.0 + s * s * s * s));
	//return -250 * exp(- 5 *s);
};

double test12Rs(double x, double y, double alpha, double beta, double sx, double sy)
{
	if (x < 1.0)
	{
		return beta;
	}
	else
	{
		return 0;
	};
};

double test12NL(double s,double x, double y, double alpha, double beta, double sx, double sy)
{
	if (x < 1.0)
	{
		return AirNL(s);
	}
	else
	{
		return IronNL(s);
	};
};

double test12Deriv(double s,double x, double y, double alpha, double beta, double sx, double sy)
{
	if (x < 1.0)
	{
		return AirDeriv(s);
	}
	else
	{
		return IronDeriv(s);
	};
};



double phcoef()
{
	return exp(1.0 - 1.0 / (physicalness * physicalness));
};

double SimpleRs(double x, double y, double alpha, double beta, double sx, double sy)
{
	if (x >= sx)
	{
		y += sy;
	};
	if (y >= 1)
	{
		y-=1.0;
	};

	double brick = 1.0/4.0;

	x /= brick;
	y /= brick;

        if (x > 1.0 && x < 2.0 && y > 1.0 && y < 3.0)
        {
                return 0;
        };

        if (x > 6.0 && x < 7.0 && y > 1.0 && y < 3.0)
	{
		return 0;
	};

	if (x > 3.0 + alpha && x < 4.0 + alpha && y > 1.0 && y < 3.0)
	{
		return beta;
	};

        return 0;
};

double SimpleNL(double s,double x, double y, double alpha, double beta, double sx, double sy)
{
	if (x >= sx)
	{
		y += sy;
	};
	if (y >= 1)
	{
		y-=1.0;
	};

	double brick = 1.0/4.0;

	x /= brick;
	y /= brick;

        if (x > 1.0 && x < 2.0 && y > 1.0 && y < 3.0)
        {
                return ArmNL(s) * phcoef() + AirNL(s) * (1.0 - phcoef());
        };

        if (x > 6.0 && x < 7.0 && y > 1.0 && y < 3.0)
	{
		return IronNL(s) * phcoef() + AirNL(s) * (1.0 - phcoef());
	};

	if (x > 3.0 + alpha && x < 4.0 + alpha && y > 1.0 && y < 3.0)
	{
		return WoodNL(s);
	};

        return AirNL(s);
};

double SimpleDeriv(double s,double x, double y, double alpha, double beta, double sx, double sy)
{

	if (x >= sx)
	{
		y += sy;
	};
	if (y >= 1)
	{
		y-=1.0;
	};

	double brick = 1.0/4.0;

	x /= brick;
	y /= brick;

        if (x > 1.0 && x < 2.0 && y > 1.0 && y < 3.0)
        {
                return ArmDeriv(s) * phcoef() + AirDeriv(s) * (1.0 - phcoef());
        };

        if (x > 6.0 && x < 7.0 && y > 1.0 && y < 3.0)
	{
		return IronDeriv(s) * phcoef() + AirDeriv(s) * (1.0 - phcoef());
	};

	if (x > 3.0 + alpha && x < 4.0 + alpha && y > 1.0 && y < 3.0)
	{
		return WoodDeriv(s);
	};

        return AirDeriv(s);
};

double MagnetRs(double x, double y, double alpha, double beta, double sx, double sy)
{
	double brick = 1.0 / 7.0;

	if (x >= sx)
	{
		y += sy;
	};
	if (y >= 1)
	{
		y-=1.0;
	};

	if (x > 6.0 * brick + alpha && x < 9.0 * brick + alpha && y > 3.0 * brick && y < 4.0 * brick)
	{
		return beta;
	};
};

double MagnetNL(double s, double x, double y, double alpha, double beta, double sx, double sy)
{

	double brick = 1.0 / 7.0;

	if (x >= sx)
	{
		y += sy;
	};
	if (y >= 1)
	{
		y-=1.0;
	};


	//armature
	if (x > 3.0 * brick && x < 4.0 * brick && y > brick && y < brick * 6.0)
	{
		return ArmNL(s);
	};

	//core, 1
	if (x > 5.0 * brick && x < 11.0 * brick && y > 1.0 * brick && y < 2.0 * brick)
	{
		return IronNL(s);
	};

	//core, 2
	if (x > 5.0 * brick && x < 11.0 * brick && y > 5.0 * brick && y < 6.0 * brick)
	{
		return IronNL(s);
	};

	//core, 3
	if (x > 10.0 * brick && x < 11.0 * brick && y > 2.0 * brick && y < 5.0 * brick)
	{
		return IronNL(s);
	};

	return AirNL(s);
};

double MagnetDeriv(double s, double x, double y, double alpha, double beta, double sx, double sy)
{
	double brick = 1.0 / 7.0;

	if (x >= sx)
	{
		y += sy;
	};
	if (y >= 1)
	{
		y-=1.0;
	};

	//armature
	if (x > 3.0 * brick && x < 4.0 * brick && y > brick && y < brick * 6.0)
	{
		return ArmDeriv(s);
	};

	//core, 1
	if (x > 5.0 * brick && x < 11.0 * brick && y > 1.0 * brick && y < 2.0 * brick)
	{
		return IronDeriv(s);
	};

	//core, 2
	if (x > 5.0 * brick && x < 11.0 * brick && y > 5.0 * brick && y < 6.0 * brick)
	{
		return IronDeriv(s);
	};

	//core, 3
	if (x > 10.0 * brick && x < 11.0 * brick && y > 2.0 * brick && y < 5.0 * brick)
	{
		return IronDeriv(s);
	};

	return AirDeriv(s);
};

double EngineNL(double s, double x, double y, double sx, double sy)
{
	double brick = 1.0 / 12.0;

	if (x >= sx)
	{
		y += sy;
	};
	if (y >= 1)
	{
		y-=1.0;
	};

	x /= brick;
	y /= brick;

	if (x < 10.0)
	{
		if ((x > 6.0 && x < 9.0 && y > 1.0 && y < 5.0) || (x > 6.0 && x < 9.0 && y > 7.0 && y < 11.0))
		{
			//magnets
			return AirNL(s);
		}
		else
		{
			//rotor
			return IronNL(s);
		};
		
	};

	if (x > 14.0)
	{
		if (x > 15.0 && x < 21.0 && y > 3.0 && y < 6.0)
		{
			//coil 1
			return AirNL(s);
		};	
		if (x > 15.0 && x < 21.0 && y > 7.0 && y < 10.0)
		{
			//coil 2 
			return AirNL(s);
		};	
		if (x > 15.0 && x < 21.0 && (y > 11.0 || y < 2.0))
		{
			//coil 3
			return AirNL(s);
		};	

		if (x < 15.0 && y > 0.0 && y < 1.0)
		{
			//air 1
			return AirNL(s);
		};

		if (x < 15.0 && y > 4.0 && y < 5.0)
		{
			//air 2
			return AirNL(s);
		};

		if (x < 15.0 && y > 8.0 && y < 9.0)
		{
			//air 3
			return AirNL(s);
		};

		//stator material
		return AirNL(s);
	};

	//air
	return AirNL(s);
};

double EngineDeriv(double s, double x, double y, double sx, double sy)
{
	double brick = 1.0 / 12.0;

	if (x >= sx)
	{
		y += sy;
	};
	if (y >= 1)
	{
		y-=1.0;
	};

	x /= brick;
	y /= brick;

	if (x < 10.0)
	{
		if ((x > 6.0 && x < 9.0 && y > 1.0 && y < 5.0) || (x > 6.0 && x < 9.0 && y > 7.0 && y < 11.0))
		{
			//magnets
			return AirDeriv(s);
		}
		else
		{
			//rotor
			return IronDeriv(s);
		};
		
	};

	if (x > 14.0)
	{
		if (x > 15.0 && x < 21.0 && y > 3.0 && y < 6.0)
		{
			//coil 1
			return AirDeriv(s);
		};	
		if (x > 15.0 && x < 21.0 && y > 7.0 && y < 10.0)
		{
			//coil 2 
			return AirDeriv(s);
		};	
		if (x > 15.0 && x < 21.0 && (y > 11.0 || y < 2.0))
		{
			//coil 3
			return AirDeriv(s);
		};	

		if (x < 15.0 && y > 0.0 && y < 1.0)
		{
			//air 1
			return AirDeriv(s);
		};

		if (x < 15.0 && y > 4.0 && y < 5.0)
		{
			//air 2
			return AirDeriv(s);
		};

		if (x < 15.0 && y > 8.0 && y < 9.0)
		{
			//air 3
			return AirDeriv(s);
		};

		//stator material
		return AirDeriv(s);
	};

	//air
	return AirDeriv(s);
};

double EngineRs(double x, double y, double beta, double sx, double sy)
{
	double brick = 1.0 / 12.0;

	if (x >= sx)
	{
		y += sy;
	};
	if (y >= 1)
	{
		y-=1.0;
	};

	x /= brick;
	y /= brick;

	if (x < 10.0)
	{
		if (x > 6.0 && x < 9.0 && y > 1.0 && y < 5.0)
		{
			//magnet 1
			return -beta * 0.5;
		};
		if (x > 6.0 && x < 9.0 && y > 7.0 && y < 11.0)
		{
			//magnet 2
			return beta * 0.5;
		};

		return 0;		
	};

	if (x > 14.0)
	{
		if (x > 15.0 && x < 21.0 && y > 3.0 && y < 6.0)
		{
			//coil 1
			return beta;
		};	
		if (x > 15.0 && x < 21.0 && y > 7.0 && y < 10.0)
		{
			//coil 2 
			return beta;
		};	
		if (x > 15.0 && x < 21.0 && (y > 11.0 || y < 2.0))
		{
			//coil 3
			return - beta;
		};	

		if (x < 15.0 && y > 0.0 && y < 1.0)
		{
			//air 1
			return 0;
		};

		if (x < 15.0 && y > 4.0 && y < 5.0)
		{
			//air 2
			return 0;
		};

		if (x < 15.0 && y > 8.0 && y < 9.0)
		{
			//air 3
			return 0;
		};

		//stator material
		return 0;
	};

	//air
	return 0;
};



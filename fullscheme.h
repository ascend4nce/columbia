
struct solverparams
{
	int n;
	double h;
	int varnum;

	int nlnum;

	double shifty, shiftx;

	double alpha, beta;

	double scaling;

	//in-region boundary parameters
	int halfnum, shiftnum;

	solverparams(int nset, double alphaset, double betaset, double shiftx, double shifty, double scaling);

	double cfunc(double s, double x, double y);

	double cfuncderiv(double s, double x, double y);
	
	double rs(double x, double y);

	double solution(double x, double y);

	double fluxx(double x, double y);

	double fluxy(double x, double y);

	int pnum(int i, int j);

	int exnum(int i, int j);

	int eynum(int i, int j);

	int p2num(int i, int j);

	int ex2num(int i, int j);

	int ey2num(int i, int j);

	int shiftindex(int);
	
	int unrollindex(int);

	void vecshiftrotate(double * values);

	void vecshiftunroll(double * values);

};
int modL(int num, int loop);
double oneresidlin(double * pt, void * params, int num);

double oneresidnl(double * pt, void * params, int num);

void fullsolverresidual(double * pt, double * ev, void * params);

matrixCSR * fullsolverderivnum(double * pt, void * params);

matrixCSR * fullsolverderivCSR(double * pt, void * params);

matrixCSR * partderivCSR(double * pt, void * params, int rank, int * ia1);


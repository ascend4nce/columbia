
void partsolverresidual(double * pt, double * ev, void * params);

matrixCSR * partsolverderivCSR(double * pt, void * params);
matrixCSR * partsolverderivnum(double * pt, void * params);

struct partsolverparams
{
	solverparams * prs;

	int rank;
	int current_rank;
	int rank2;
	int rank9;

	int interpdata;

	int localshiftnum;
	
	int * ia9;
	int * ia1;

	double * VTAV;
	double * Aepsex;
	double * Aepsv;
	double * VTUSU;
	double * I9V;

	double * alphaexample;
	double * betaexample;
	double * shiftexample;
	double * interpmatrix;

	bool isLoaded;
	bool isInterpLoaded;
	double * V;
	double * U;

	void testnl(int rk, int ntests);

	void findoptimal(int rk, double * vec, double * optapprox);

	partsolverparams(int rankset, int rank2set, char * basis, solverparams * prsset);
	
	void applyreductionshift(int);
	
	void uploadbasis(char *);

	void uploadinterpbasis(char *);

	~partsolverparams();
};

void reducedtest(partsolverparams * PPRS, int rk, double * smallpoint);
void reducedtestmark2(partsolverparams * PPRS, int rk, double * smallpoint);
void restoresolution(partsolverparams * PPRS, double * smallpoint, double * point);

class matrixCSR
{
public:
	int m,n;
	double * a;
	int * ia;
	int * ja;

	matrixCSR();

	~matrixCSR();

	void print();
	
	matrixCSR(int mset, int nset, double * aset, int * iaset, int * jaset);

	void matvec(double * vec, double * res);
	
	void leftmatvec(double * vec, double * res);
	
	void solveprecgmres(double * rs, double * x, double lstolerance);
 
	void solvepreccovcg(double * rs, double * x, double lambda, double lstolerance);
 
	void solvesvd(double * rs, double * x);

	void solvecovsvd(double * res, double * x, double lambda);
 

};

#define METHOD_PRECGMRES 907
#define METHOD_SVD 908

#define METHOD_PRECCOVCG 1017
#define METHOD_COVSVD 1019

typedef void (* eval)(double *, double *, void *);
typedef matrixCSR * (*derivCSR)(double *, void *);

void newton(int n, double * point, eval function, derivCSR derivmaker, void * params, double tolerance, int maxit, double startsigma, double lstolerance, int method, int covmethod);
void newtonmark2(int n, double * point, eval function, derivCSR derivmaker, void * params, double tolerance, int maxit, double startsigma, double lstolerance, int method);


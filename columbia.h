
void simpletest(solverparams * , double *);

void builddata(int gridN, int as, int bs, int ss, int rankmax, int basenum, int snapshots, int nlsnapshots, char * basis, char * assist, char * interp);

void printreference(solverparams * , double *, double * , char * );

void putsolution(solverparams * , double * );

void printrs(solverparams * , char * );

void printsimple(solverparams *, double *, char *, char *, char *);

void printfluxmodule(solverparams *, double *, char *, char *);

void multitest(solverparams * PRSinit, double *, int doubling);

solverparams * gridtoparams(int gridN,int a, int b, int s, int as, int bs, int ss);

void printdependencynodes(partsolverparams * PPRS, char * where1, char * where9);


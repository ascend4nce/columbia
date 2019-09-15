#define printlevel 1

//method constants
#define fullnewtonmethod 2 //'1' is the trust-region enhanced Newton methods with optimization parameter 'nu'; should be more accurate than '2' in exact arithmetics, but uses CG for covariation matrices, suffers from matrix being ill-conditioned, so should not be used. '2' is the broken line approximate method that only uses GMRES for jacobian matrix, usually dramatically outperforms '1'.
#define fullderivatives 1 //'1' means analytic expressions are used for derivatives and derivatives are stored in the sparse matrix CSR format. '2' means numerical derivatives are used; was used for testing only; they are approximate, slow and memory-consuming (because a lot of zeros are stored.); for '2', GRIDsize should be <= 40, otherwise memory errors are likely due to memory limit.
#define reducednewtonmethod 2 //same as above, but for reduced smaller problem
#define reducedderivatives 1 //same as above, but for reduced smaller problem

//gmres and cg constants
#define maxgmresit 50 //GMRES iterations limit (used in Newton-2)
#define maxcgit 100 //CG iterations limit (used in Newton-1)
#define gmrestart 50 //GMRES iterations before first restart
//it seems that preconditioner works fine for GMRES, so restarts are not needed

//newton constants
#define Tolerance 1.0e-10 //Newton iterations tolerance
#define LStolerance 1.0e-12 //GMRES/CG (that is used inside Newton iterations) tolerance 
#define MAXnewtoniterations 1000 //Newton iteration limit; around 10 iterations are enough if trust regions are not used and around 50 should be enough if they are used.

//grid size
#define GRIDsize 40 //'N' in the reports; defines vertical grid size, horizontal is then equal to 2N
#define GridNumber 2 //Number of hierarhical grids, that are used for obtaining better Newton starting point. GridNumber >= 1; if GridNumber = k, then GRIDsize should be divisible by 2^(k-1), otherwise errors will occur

//scaling for linear part. Set around 1.0e-6 for physical values
#define Scaling 1.0e-6 //Scaling factor that is applied to equations with nontrivial nonlinearity. Should not be changed if realistic material properties are used, works good.

#define physicalness 1.0 //A constant that smoothly pushes the nonlinearity from physical values to a linear function as it falls from 1 to 0. Should be equal to 1.0; was used for testing only
#define Shiftx (6.0/10.0) //X axis position of the shift vertical line. Should be aligned with grid square edges: value of 7.0/15.0 is illegal if GRIDsize = 10.Range: [0;2); changing it may break geometry (rip off some rectangles).

//test constants
#define Alpha (0.0 / 12.0) //horizontal coil shift. Range [0; 1). Not used in final experiments due to bad ranks (see stage 2 report)
#define Beta 9.0e+5 //Right hand side values. Nontrivial values are ranged somewhere in [1.0e+4, 5.0e+6).
#define Shifty (3.0 / 20.0) //Vertical cyclic shift value. Range [0;1). Should be aligned with grid square edges: value of 7.0/15.0 is illegal if GRIDsize = 10.

//basis building parameters
#define alphagrid 1 //alpha parameter uniform splitting. Should normally be 1 cause the parameter alpha is 'bad' (see stage 2 report for details)
#define betagrid 20 //beta parameter uniform splitting in logarithmic scale over interval [3.0e+4, 3.0e+6). see columbia.cpp -> gridtoparams for details
#define shiftgrid 20 //cyclic shift parameter uniform splitting over interval [0,1); GRIDsize should be divisible by shiftgrid
#define RANKnonlinearity 300 //number of interpolation nodes == dimensionality of nonlinearity approximation subspace. Should be >> RANKreduction, also depends on the GRIDsize.
#define RANKmaxreduction 50 //maximum reduction rank value: for this value, corresponding parts of auxiliary matrices will be computed and stored. Increasing it with a fixed RANKnonlinearity value will reduce the nonlinearity approximation performance. Ranks above 50 usually are useless, normally used values are around 20-45. Ranks being too high may result in Newton iterations not converging. 
#define Snapshots 80 //number of snapshots used for constructing the approximation basis. Snapshots > RANKmaxreduction.
#define NLSnapshots 1000 //number of nonlinear function basis snapshots; is normally const * RANKnonlinearity; must be > RANKnonlinearity; too high values may result in memory problems, cause basis building will require an SVD of a matrix of size NLSnapshots * GRIDsize in the offline stage


//reduction constants
#define RANKreduction 40 //reduction rank for the current test

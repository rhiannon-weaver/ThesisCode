#ifndef POISSONEM
#define POISSONEM
#include <math.h>
#include <vector>
#include <gsl/gsl_randist.h>
#define EPSILON 0.0000000000000000001
extern std::ostream cnull;
using namespace std;
double LogMixturePoissonLike(vector<int>*x, vector<double>*lambda, vector<double>*pi, int return_BIC = 1);
int PoissonEMModel(vector<int>*x, vector<double>*lambda_start, vector<double>*pi_start, double* objective, int max_iter=5000);
int DataDrivenStartValRefine(vector<int>* flows, vector<short int>* state, double* q, double* alpha, double* omega, double* rho, double* gamma, double* nu, ostream& printvalues=cnull);


#endif

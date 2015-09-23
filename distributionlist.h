#ifndef DISTRIBUTIONLIST
#define DISTRIBUTIONLIST
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <time.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <cstdarg>
#include <string>
#define NUM_DISTRIBUTIONS 10

using namespace std;

/*=========================

0: normal(mean, stdev)
1: logistic
2: beta       
3: uniform(min, max)
4: gamma(a, b) mu=ab
5: weibull
6: truncated normal(mu, sd, min, max)
7: scaled beta(min, max, alpha, beta)
8: multinomial
9: finitediscrete
10: poisson(lambda)
==========================*/

class DistributionList{
  
 public:
  int GetDistributionIntVal(const char* name);
  DistributionList(int seedval = 0);
  void Simulate(int which, int num, ...);
  double LogDensity(int which, ...);
  void Seed(int seedval);
  ~DistributionList();
  
 private:
 
  vector<string> names;
  time_t T;
  gsl_rng* random_number_generator; 
  const gsl_rng_type** generator_types;
};
#endif

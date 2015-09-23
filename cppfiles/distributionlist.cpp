#include "../distributionlist.h"
#define GAMMA_MIN 0.00000000001
DistributionList::DistributionList(int seedval)
{
  generator_types = gsl_rng_types_setup();
  random_number_generator = gsl_rng_alloc(gsl_rng_mt19937);
  if(seedval == 0)
    gsl_rng_set(random_number_generator, time(&T));
  else
    gsl_rng_set(random_number_generator, seedval);

  names.push_back("normal");     //0
  names.push_back("logistic");   //1
  names.push_back("beta");       //2
  names.push_back("uniform");    //3
  names.push_back("gamma");      //4
  names.push_back("weibull");    //5
  names.push_back("truncnorm");  //6
  names.push_back("scaledbeta"); //7
  names.push_back("multinomial");//8
  names.push_back("finitediscrete");//9
  names.push_back("poisson"); //10

}

int DistributionList::GetDistributionIntVal(const char* name)
{ 
  int i, found;
  i = 0;
  found = 0;
  while(!found && i < NUM_DISTRIBUTIONS)
    {
      if( strcmp(name, names[i].c_str()) == 0 )
        found = 1;
      else
        i++;
    }
  if(!found)
    {
      i = -1;
      cout << "Distribution " <<  name << " not found" << endl;
    }
  return i;
  
}


void DistributionList::Simulate(int which, int num, ...)
{
  int i;
  va_list ap;
  double* vals_dbl_ptr;
  unsigned int* vals_int_ptr;
  vector<double> par;
  gsl_ran_discrete_t * ran_preproc;
  double* pk;
  va_start(ap, num);
  switch(which)
    {
    case 0:
      //normal
      vals_dbl_ptr = va_arg(ap, double*);
      par.push_back(va_arg(ap, double));
      par.push_back(va_arg(ap, double));
      for(i = 0; i < num; i++)
        vals_dbl_ptr[i] = gsl_ran_gaussian(random_number_generator, par[1]) + par[0];
      break;
    case 1:
      //logistic  
      vals_dbl_ptr = va_arg(ap, double*);
      par.push_back(va_arg(ap, double));
      par.push_back(va_arg(ap, double));
      for(i = 0; i < num; i++)
        vals_dbl_ptr[i] = gsl_ran_logistic(random_number_generator, par[1]) + par[0];
      break;
    case 2:
      //beta
      vals_dbl_ptr = va_arg(ap, double*);
      par.push_back(va_arg(ap, double));
      par.push_back(va_arg(ap, double));
      for(i = 0; i < num; i++)
        vals_dbl_ptr[i] = gsl_ran_beta(random_number_generator, par[0], par[1]);
      break;
    case 3:
      //uniform
      vals_dbl_ptr = va_arg(ap, double*);
      par.push_back(va_arg(ap, double));
      par.push_back(va_arg(ap, double));
      for(i = 0; i < num; i++)
        vals_dbl_ptr[i] = gsl_ran_flat(random_number_generator, par[0], par[1]);
      break;
    case 4:
      //gamma
      vals_dbl_ptr = va_arg(ap, double*);
      par.push_back(va_arg(ap, double));
      par.push_back(va_arg(ap, double));
      for(i = 0; i < num; i++)
        vals_dbl_ptr[i] = gsl_ran_gamma(random_number_generator, par[0], par[1]);
      break;
    case 5:
      //weibull
      vals_dbl_ptr = va_arg(ap, double*);
      par.push_back(va_arg(ap, double));
      par.push_back(va_arg(ap, double));
      for(i = 0; i < num; i++)
        vals_dbl_ptr[i] = gsl_ran_weibull(random_number_generator, par[0], par[1]);
      break;
    case 6:
      //truncated normal
      vals_dbl_ptr = va_arg(ap, double*);
      par.push_back(va_arg(ap, double));
      par.push_back(va_arg(ap, double));
      par.push_back(va_arg(ap, double));
      par.push_back(va_arg(ap, double));
      for(i = 0; i < num; i++)
	{
	  vals_dbl_ptr[i]= gsl_ran_gaussian(random_number_generator, par[1]) + par[0];
	  while( (vals_dbl_ptr[i] < par[2]) || (vals_dbl_ptr[i] > par[3]))
	    vals_dbl_ptr[i]= gsl_ran_gaussian(random_number_generator, par[1]) + par[0];
	}
      break;
    case 7:
      //scaled beta
      vals_dbl_ptr = va_arg(ap, double*);
      par.push_back(va_arg(ap, double)); //minimum bound
      par.push_back(va_arg(ap, double)); //maximum bound
      par.push_back(va_arg(ap, double)); //alpha
      par.push_back(va_arg(ap, double)); //beta
      for(i = 0; i < num; i++)
        vals_dbl_ptr[i] = (par[1]- par[0])*gsl_ran_beta(random_number_generator, par[2], par[3])+par[0];
      break;
    case 8:
      //multinomial, assuming end-to-end in vector
      vals_int_ptr = va_arg(ap, unsigned int*);
      par.push_back(va_arg(ap, double)); //N
      par.push_back(va_arg(ap, double)); //K
      pk = va_arg(ap, double*);
      gsl_ran_multinomial(random_number_generator, (size_t)(par[1]),(unsigned int)(par[0]), pk, vals_int_ptr);
      break;
    case 9:
      //general discrete
      vals_int_ptr = va_arg(ap, unsigned int*);
      ran_preproc = va_arg(ap, gsl_ran_discrete_t*); //preprocessed lookup table
      for(i = 0; i < num; i++)
	vals_int_ptr[i] = gsl_ran_discrete(random_number_generator, ran_preproc);
      break;
    case 10:
      //poisson
      vals_dbl_ptr = va_arg(ap, double*);
      par.push_back(va_arg(ap, double));
      for(i = 0; i < num; i++)
        vals_dbl_ptr[i] = gsl_ran_poisson(random_number_generator, par[1]);
      break;
    }
  va_end(ap);
  
}


double DistributionList::LogDensity(int which, ...)
{
  double val;
  double x;
  unsigned int x_int;
  unsigned int * x_int_ptr;
  vector<double> par;
  double* x_dbl_ptr;
  double* pk;
  va_list ap;
  va_start(ap, which);
  switch(which)
    {
      //normal: par2 is standard deviation
      //gamma: par2 is the denominator value (scale, not rate)
    case 0:
      //normal
      x = va_arg(ap, double);
      par.push_back(va_arg(ap, double));
      par.push_back(va_arg(ap, double));
      val = gsl_ran_gaussian_pdf(x - par[0], par[1]);
      break;
    case 1:
      //logistic
      x = va_arg(ap, double);
      par.push_back(va_arg(ap, double));
      par.push_back(va_arg(ap, double));
      val = gsl_ran_logistic_pdf(x - par[0], par[1]); 
      break;
    case 2:
      //beta
      x = va_arg(ap, double);
      par.push_back(va_arg(ap, double));
      par.push_back(va_arg(ap, double));
       val = gsl_ran_beta_pdf(x, par[0], par[1]); 
      break;
    case 3:
      //uniform
      x = va_arg(ap, double);
      par.push_back(va_arg(ap, double));
      par.push_back(va_arg(ap, double));
      val = gsl_ran_flat_pdf(x, par[0], par[1]); 
      break;
    case 4:
      x = va_arg(ap, double);
      par.push_back(va_arg(ap, double));
      par.push_back(va_arg(ap, double));
      if(x < 0.0)
	{
	 x = GAMMA_MIN;
        }
      val = gsl_ran_gamma_pdf(x, par[0], par[1]);
      //gamma
      break;
    case 5:
      x = va_arg(ap, double);
      par.push_back(va_arg(ap, double));
      par.push_back(va_arg(ap, double));
      //weibull
      val = gsl_ran_weibull_pdf(x, par[0], par[1]);
      break;
    case 6:
      x = va_arg(ap, double);
      par.push_back(va_arg(ap, double));
      par.push_back(va_arg(ap, double));
      par.push_back(va_arg(ap, double));
      par.push_back(va_arg(ap, double));
      par.push_back(gsl_cdf_gaussian_P(par[3]-par[0], par[1]) - gsl_cdf_gaussian_P(par[2]-par[0],par[1]));
      val = (((x >= par[2]) && (x <= par[3]))*(gsl_ran_gaussian_pdf(x - par[0],par[1])))/par[4];
      break;
    case 7:
      x = va_arg(ap, double);
      par.push_back(va_arg(ap, double)); //lower bound
      par.push_back(va_arg(ap, double)); //upper bound
      par.push_back(va_arg(ap, double)); //alpha
      par.push_back(va_arg(ap, double)); //beta
      val = (1/(par[1] - par[0] ) )*gsl_ran_beta_pdf( (x-par[0])/(par[1] - par[0] ), par[2], par[3]) * ((x>=par[0]) && (x <= (par[1])));
      break;
    case 8:
      x_int_ptr = va_arg(ap, unsigned int*);
      par.push_back(va_arg(ap, double)); //K
      pk = va_arg(ap, double*);
      val = gsl_ran_multinomial_pdf(par[0], pk, x_int_ptr);
      break;
    case 9:
      //general discrete
      x_int = va_arg(ap, unsigned int);
      x_dbl_ptr = va_arg(ap, double*);
      val = par[x];
      break;
    case 10:
      //poisson
      x_int = va_arg(ap, int);
      par.push_back(va_arg(ap, double));
      val = gsl_ran_poisson_pdf(x, par[1]);
      break;
    }
  va_end(ap);
  return(log(val));
}

DistributionList::~DistributionList()
{
  int i;
  gsl_rng_free(random_number_generator);
 
}

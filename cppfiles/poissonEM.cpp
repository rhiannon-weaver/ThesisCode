#include "../poissonEM.h"
#include "../datparobjects.h"
std::ostream cnull(0);

double LogMixturePoissonLike(vector<int>*x, vector<double>*lambda, vector<double>* pi, int return_BIC)
{
  double loglike = 0.0;
  double sumprod, cur_lam, cur_pi;
  int i, j, dim, n;
  int cur_x;
  n = x->size();
  dim = lambda->size();
  for(i = 0; i < n; i++)
    {
      cur_x = (int)x->operator[](i);
      sumprod = 0;
      for(j = 0; j < dim; j++)
	{
	  cur_lam = lambda->operator[](j);
	  cur_pi = pi->operator[](j);
	  sumprod += cur_pi*gsl_ran_poisson_pdf(cur_x, cur_lam); 
	}
      loglike += log(sumprod);
    }
  if(return_BIC)
    loglike = -2.0*loglike + (2.0*(double)(lambda->size()) - 1.0)*log((double)(n));
  return(loglike);
}

int PoissonEMModel(vector<int>* x, vector<double>* lambda_start, vector<double>* pi_start, double* objective,  int max_iter)
{
  int i, j;

  int numgroups = lambda_start->size();
  int numobs = x->size();
  int converged = 0;
  int steps = 0;

  int cur_x;

  //length = numgroups
  vector<double> ownership;
  double denom;

  vector<double> pri_lambda;
  vector<double> pri_pi;
  double pri_objective = LogMixturePoissonLike(x, lambda_start, pi_start, 0);

  vector<double> new_lambda;
  vector<double> new_pi;
  double new_objective;

  for(j = 0; j < numgroups; j++)
    {
      ownership.push_back(0.0);
      pri_lambda.push_back(lambda_start->operator[](j));
      pri_pi.push_back(pi_start->operator[](j));

      new_lambda.push_back(0.0);
      new_pi.push_back(0.0);
    }

  while( (!converged) && (steps < max_iter))
    {
      for(j = 0; j < numgroups; j++)
	{
	  new_pi[j] = 0.0;
	  new_lambda[j] = 0.0;
	}
      for(i = 0; i < numobs; i++)
	{
	  cur_x = x->operator[](i);
	  denom = 0.0;
	  for(j = 0; j < numgroups; j++)
	    {
	      ownership[j] = pri_pi[j]*gsl_ran_poisson_pdf(cur_x, pri_lambda[j]);
	      denom += ownership[j];
	    }
	  for(j = 0; j < numgroups; j++)
	    ownership[j] /= denom;
	  for(j = 0; j < numgroups; j++)
	    {
	      new_pi[j] += ownership[j];
	      new_lambda[j] += (double)cur_x*ownership[j];
	    }
	}
      //now we have unnormalized bits for lambda, pi
      denom = 0;
      for(j = 0; j < numgroups; j++)
	{
	  new_lambda[j] /= new_pi[j];
	  denom += new_pi[j];
	}
      for(j = 0; j < numgroups; j++)
	new_pi[j] /= denom;

      new_objective = LogMixturePoissonLike(x, &new_lambda, &new_pi, 0);
      if(fabs(new_objective - pri_objective) < EPSILON)
	{
	  converged = 1;
	}
      else
	{
	  steps++;
	  pri_objective = new_objective;
	}
    } 
  if(!converged)
    steps = max_iter +1;

  *objective = new_objective;
  for(j =0; j < numgroups; j++)
    {
      lambda_start->operator[](j) = new_lambda[j];
      pi_start->operator[](j) = new_pi[j];
    }
  return(steps);
}



int DataDrivenStartValRefine(vector<int>* flows, vector<short int>* state, double* q, double* alpha, double* omega, double* rho, double* gamma, double* nu, ostream& printvalues)
{
  //Assuming that the off states have already been labeled, and off_lambda has already been set, ie
  //that InternalSet and SetInitialStartState have been done already.

  int i, converged;
  int num_init_spike, num_init_decay;
  int sum_spike_vals, sum_decay_vals;
  int init_spike_sum, init_decay_sum;
  
  vector<int> x;
  vector<double> lambda_start_1par;
  vector<double> pi_start_1par;

  vector<double> lambda_start_2par;
  vector<double> pi_start_2par;
  
  double BIC_1par;
  double BIC_2par;
  double outlier_prob;
  int max_count;
  int spike_ephemeral;

  double ownership;


  //to avoid infinities and Nan's, initialize to 1 of each, seeing a 1.
  num_init_spike = 1;
  num_init_decay = 1;
  sum_spike_vals = 1;
  sum_decay_vals = 1;

  

  for(i = 0; i < flows->size(); i++)
    {
      if(state->operator[](i) != 0)
	{
	  x.push_back(flows->operator[](i));
	  if(state->operator[](i) == 1)
	    {
	      num_init_spike++;
	      sum_spike_vals += flows->operator[](i);
	    }
	  else
	    {
	      num_init_decay++;
	      sum_decay_vals += flows->operator[](i);
	    }
	}

    }

  

  lambda_start_1par.push_back((double)(sum_spike_vals + sum_decay_vals)/(double)(x.size()+2));
  pi_start_1par.push_back(1.0);
  BIC_1par = LogMixturePoissonLike(&x, &lambda_start_1par, &pi_start_1par, 1);
  
  lambda_start_2par.push_back( (double)(sum_spike_vals)/(double)(num_init_spike));
  lambda_start_2par.push_back( (double)(sum_decay_vals)/(double)(num_init_decay));

  pi_start_2par.push_back( (double)(num_init_spike)/((double)(x.size()+1)));
  pi_start_2par.push_back( 1.0 - ((double)(num_init_spike)/((double)(x.size()+1))));

  converged = PoissonEMModel(&x, &lambda_start_2par, &pi_start_2par, &BIC_2par);
  BIC_2par = -2.0*BIC_2par + (3 * log((double)x.size()));

  if(lambda_start_2par[0] > lambda_start_2par[1])
    {
      //switch them around
      outlier_prob = lambda_start_2par[0];
      lambda_start_2par[0] = lambda_start_2par[1];
      lambda_start_2par[1] = outlier_prob;
      outlier_prob = pi_start_2par[0];
      pi_start_2par[0] = pi_start_2par[1];
      pi_start_2par[1] = outlier_prob;
    }
  
  if(!(printvalues == NULL))
    {
      printvalues << x.size()-2 << " " << num_init_spike-1 << " " << num_init_decay-1 << " " <<  lambda_start_1par[0] << " " << BIC_1par << " " << lambda_start_2par[0] << " " << lambda_start_2par[1] << " " << pi_start_2par[0] << " " << pi_start_2par[1] << " " << BIC_2par << " " << converged << endl;
    }

  //First check to see that the "spike" state in the 2-par model is the rare state
  spike_ephemeral = 1;
  if(lambda_start_2par[0] > lambda_start_2par[1])
    {
      if(pi_start_2par[0] > pi_start_2par[1])
	spike_ephemeral = 0;
    }
  else
    {
      if(pi_start_2par[0] < pi_start_2par[1])
	spike_ephemeral = 0;
    }

  if( (BIC_2par < BIC_1par) && (spike_ephemeral))
    {
      //Use the 2-parameter settings to initialize the MCMC
      *q = lambda_start_2par[0];
      *omega = lambda_start_2par[1]/lambda_start_2par[0];
      *alpha = 0.65;  //middling, IDEK
      //okay now we have low prob = 0th position, high prob = 1st position
      //if pr(zi = spike | xi) > 0.75 we'll set it as a spike (to account for the alpha)
      for(i = 0; i < flows->size(); i++)
	{
	  if(state->operator[](i) != 0)
	    {
	      ownership = (pi_start_2par[1]*gsl_ran_poisson_pdf((int)flows->operator[](i), lambda_start_2par[1]))/( (pi_start_2par[1]*gsl_ran_poisson_pdf((int)flows->operator[](i), lambda_start_2par[1]))+ (pi_start_2par[0]*gsl_ran_poisson_pdf((int)flows->operator[](i), lambda_start_2par[0])));
	      if(ownership > 0.75)
		state->operator[](i) = 1;
	      else
		state->operator[](i) = 2;
	    }
	}
      //setting rho/gamma/nu

    }

  else
    {
      //Use the 1-parameter settings to initialize the MCMC
      //Re-parameterize the psi parameterization for decay states to do decay-to-spike gamma 
      //small, and time dependent between decay-off and decay-decay

      //IF WE DO THIS WE HAVE TO OUTPUT THE PSI PARAMETERIZATION FOR EACH MACHINE!
      //This isn't too difficult to do though.
      //But it would be easier to just set it to something different.  Let's do that.

      //use the initial 1-parameter lambda estimate for q
      *q = lambda_start_1par[0];
      //Set alpha as low, indicating that a spike, if it occurs, is ephemeral
      *alpha = 0.25;
      //Set omega large (there may be spikes in this profile but we aren't seeing them)
      *omega = 6.0;

      //Set the spike to anything not spike parameter high, anything->spike low.
      gamma[1] = 0.001;
      //Set the others
      gamma[0] = 0.1; //doesn't really matter
      gamma[2] = (double)(x.size())/(double)(flows->size());  //should be ratio of off to decay
      rho[2] = 0.01;
      nu[2] = 15.0;
      
      rho[0] = 0.1;
      nu[0] = 15.0;

      rho[1] = 0.8;
      nu[1] = 5.0;
      //also need to change the priors?  (ans: no, they are uniform)
      
      max_count = -1;
      outlier_prob = 0.0;
      while(outlier_prob < 0.995)
	{
	  max_count++;
	  outlier_prob += gsl_ran_poisson_pdf(max_count, lambda_start_1par[0]);
	}
      for(i = 0; i < flows->size(); i++)
	{
	  if(state->operator[](i) != 0)
	    {
	      if(flows->operator[](i) > (short int) max_count)
		state->operator[](i) = 1;
	      else
		state->operator[](i) = 2;
	    }
	}
    }
  return(converged);
}

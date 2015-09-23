#include "../helperfunctions.h"

double trunc_norm_logjump(double val1,
			  double val2,
			  double lim1, 
			  double lim2, 
			  double prop_sd, 
			  DistributionList* distributions)
{
  double x = val1;
  double curr = val2;
  return( log(((x >= lim1) & (x <= lim2))) + (distributions->LogDensity(0, x, curr, prop_sd)));
}

double trunc_norm_proposal(double current, 
			   double proposalsd, 
			   double lowerlim, 
			   double upperlim, 
			   DistributionList* distributions)
{
  double val;
  distributions->Simulate(6, 1, &val, current, proposalsd, lowerlim, upperlim);
  return val;
}

double BetaConstrain(double val)
{
  double newval;
  newval = val;
  if(newval > BETA_MAX)
    newval = BETA_MAX;
  if(newval < BETA_MIN)
    newval = BETA_MIN;
  return newval;
}

double sinecurveprob_single(short int pristate, 
			    short int nextstate,
			    short int t,
			    double rho,
			    double gamma,
			    double nu, 
			    curr_par_obj* current_values)
{

  double* prob_mult;
  double prob1;
  double gamprob;

  //adaptation for consecutive off states: just return 1.0
  if( (pristate == 0) && (nextstate == 0))
    {
      prob1 = 1.0;
      gamprob = 1.0;
    }
  else //we do not have a 0->0 transition
    {
      prob1 = (rho/(nu + 2.0))*(sin(2.0*M_PI*t/24.0) + nu + 1.0);
      prob_mult = current_values->psi_parameterization[pristate][nextstate];
      //cout << endl << "sinecurveprob_single: prob1 = " << prob1 << ", gamma = " << gamma << " prob_mult == " << prob_mult[0] << " " << prob_mult[1] << endl; 
      if((prob_mult[0] == 0) && (prob_mult[1] == 0))
	{
	  prob1 = 0.0;
	  gamprob = 0.0;
	}
      else
	{
	  if(prob_mult[0] == -1)
	    prob1 = 1.0 - prob1;
	  else 
	    prob1 = pow(prob1, prob_mult[0]);
	  
	  if(prob_mult[1] == -1)
	    gamprob = 1 - gamma;
	  else
	    gamprob = pow(gamma, prob_mult[1]);
	}
    }
  return(prob1*gamprob);
}


double sum_alpha_k(short int beginval, short int endval, curr_par_obj* current_values, double alpha, double omega)
{
  double lagsum;
  lagsum = 0.0;
  int i, end_trunc;
  if(endval >= current_values->T)
    end_trunc = endval -1;
  else
    end_trunc = endval;
  for(i = beginval; i <= end_trunc; i++)
    {
      if(current_values->Lagstate(i) != -1.0)
	lagsum += pow(alpha, current_values->Lagstate(i));
    }
  if(omega > 0.0)
    lagsum *= (omega - 1);
  return lagsum;
}


double log_prod_y_alpha(mach_dat_obj* data, curr_par_obj* current_values, double alpha, double omega, short int birthdate, short int deathdate)
{
  double logprod;
  int i, end_trunc;
  logprod = 0.0;
  if(current_values->Deathdate() >= current_values->T)
    end_trunc = current_values->T -1;
  else
    end_trunc = current_values->Deathdate();
  for(i = current_values->Birthdate(); i <= end_trunc; i++)
    {
      if(current_values->Lagstate(i) != -1)
	logprod += data->GetTotalCount(i) * log(1.0 + pow(alpha, current_values->Lagstate(i))*(omega - 1.0));
    }
  return logprod;
}

/* this has been adapted for the model that counts poisson zero states */
double NextStateProb(int pristate, int which_parameter, double value, curr_par_obj* current_values, global_par_obj* global_pars)
{
  double logpost;
  int i, end_trunc;
  double cur_rho;
  double cur_gam;
  double cur_nu;
  
  logpost = 0.0;

  switch(which_parameter)
    {
    case 0:
      cur_rho = value;
      cur_gam = current_values->gamma[pristate];
      cur_nu = current_values->nu[pristate];
      break;
    case 1:
      cur_rho = current_values->rho[pristate];
      cur_gam = value;
      cur_nu = current_values->nu[pristate];
      break;
    case 2:
      cur_rho = current_values->rho[pristate];
      cur_gam = current_values->gamma[pristate];
      cur_nu = value;
      break;
    }
  
  if(current_values->Deathdate() >= current_values->T)
    end_trunc = current_values->T - 1;
  else
    end_trunc = current_values->Deathdate();

  for(i = (current_values->Birthdate() +1); i <= end_trunc; i++)
    {

      if(!( (pristate == 0) && (current_values->State(i) == 0)))//adaptation to remove off->off transitions
	{
	  if(i != 0)
	    {
	      if(current_values->State(i-1) == pristate)
		{
		  logpost += log(sinecurveprob_single(pristate, 
						      current_values->State(i),
						      i - current_values->tstar,
						      cur_rho,
						      cur_gam,
						      cur_nu,
						      current_values));
		  if(!(current_values->Is_immune()) &&(pristate != 0))
		    {
		      if( i != current_values->Deathdate())
			logpost += log(global_pars->survival_rate[i]);
		      else
			logpost += log(1.0 - global_pars->survival_rate[i]);
		    }
		}
	    }
	}
    }
  return(logpost);
}


void PathToChange(int index, curr_par_obj* current_values, mach_dat_obj* data, int* beginend)
{
 
  int i, found, end_trunc;
  int num_off_off = 0;
  int num_off_seq = 0;
  int in_off_seq;
  found = 0;
  if(current_values->Deathdate() >= current_values->T)
    end_trunc= current_values->T -1;
  else
    end_trunc = current_values->Deathdate();
  //find end
  if(index >= end_trunc)
    {
      beginend[1] = index;
    }
  else
    {
      i = index + 1;
      if(index >= current_values->T)
	{
	  cout << "In PathToChange, index = " << index << endl;
	}
      if( (data->GetTotalCount(index) == 0) && (current_values->State(i) == 0))
	{
	  /*go until we find a spike or a decay*/
	  while( (!found) && (i < current_values->Deathdate()))
	    {
	      if(current_values->State(i) != 0)
		{
		  found = 1;
		}
	      else
		{
		  i++;
		}
	    }
	}
      else
	{
	  /*If the count is non-zero then the number of consecutive off states afterward cannot change */
	  /*Also if the count is zero but the following state is not off, then we're at the old version too*/ 
	  while( (!found) && (i <= end_trunc))
	    {
	      if(current_values->State(i) != 2)
		{
		  found = 1;
		}
	      else
		{
		  i++;
		}
	    }
	}
      beginend[1] = i;
    }

  //now for the beginning
  found = 0;
  if( index == current_values->Birthdate())
    {
      beginend[0] = current_values->Birthdate();
    }
  else //we are not at the very beginning
    {
      i = index -1;
      // cout << "In PathToChange (looking backward), index = " << index << endl;
      if( !((data->GetTotalCount(index) == 0) && (current_values->State(i) == 0)))
	{				       
	  beginend[0] = i;
	}
      else
	{
	  //go backward until we find a non-off state
	  /*go until we find a spike or a decay*/
	  while( (!found) && (i > current_values->Birthdate()))
	    {
	      if(current_values->State(i) != 0)
		{
		  found = 1;
		}
	      else
		{
		  i--;
		}
	    }
	}

      beginend[0] = i;
    }



  in_off_seq = 0;
  if(current_values->State(beginend[0]) == 0)
    {
      num_off_seq++;
      in_off_seq = 1;
    }    
  for(i = beginend[0]+1; i< beginend[1]; i++)
    {
      if( (current_values->State(i-1) == 0) && (current_values->State(i) == 0))
	num_off_off++;
      if( (current_values->State(i-1) == 0) && (current_values->State(i) != 0))
	in_off_seq = 0;
      if( (current_values->State(i-1) != 0) && (current_values->State(i) == 0))
        {
	  num_off_seq++;
	  in_off_seq = 1;
	}
    }
  beginend[2] = num_off_off;
  beginend[3] = num_off_seq;
  return;
}




double GetLogPoissonProbabilities(mach_dat_obj* data, curr_par_obj* current_values, int startval, int endval, vector<double>* indiv_pois)
{
  int i, end_trunc;
  double probability_sum;
  double temp_inf_check;
  double ratecheck;
  int not_infinity;
  double mu;

  if(endval >= current_values->Deathdate())
    {
      end_trunc = current_values->Deathdate();
      if(end_trunc >= current_values->T)
	end_trunc = current_values->T - 1;
    }
  else
    end_trunc = endval;

  probability_sum = 0;
  not_infinity = 1;
  if(startval >= current_values->Birthdate())
    i = startval;
  else
    i = current_values->Birthdate();
  if(PRINT_PATH_STEP_PROB)
    {
      while(indiv_pois->size() > 0)
	indiv_pois->pop_back();
    }
  while( (i <= end_trunc) )
    {
      if(current_values->State(i) == 0)
	{
	  if(i >= current_values->T)
	    {
	      cout << "In GetLogPoissonProbabilities, index i = " << i << endl;
	    }
	  if(data->GetTotalCount(i) == 0)
	    temp_inf_check = 0;
	  else
	    temp_inf_check = (-1.0*INFINITY);
	}
      else
	{
	  mu =  current_values->q*(1.0 + pow(current_values->alpha,current_values->Lagstate(i))*(current_values->omega - 1.0));
	  ratecheck = current_values->GetRate(i);
	  if(ratecheck != mu)
	    {
	      cout << endl<< "Rates not equal : Helperfunctions mu = " << mu << " GetRate = " << ratecheck << endl;
	    }

	  if(i >= current_values->T)
	    {
	      cout << "In GetLogPoissonProbabilities, index i = " << i << endl;
	    }
	  temp_inf_check = log(gsl_ran_poisson_pdf(data->GetTotalCount(i), mu));
          //cout << endl << "poissonprobs: mu = " << mu << ", dat = " << data->fl[i]  << ", temp_check = "<< temp_inf_check; 
	}
      if(temp_inf_check == (-1.0*INFINITY))
	{
	  not_infinity = 0;
	  probability_sum = (-1.0 * INFINITY);
	}
      if(not_infinity)
	probability_sum += temp_inf_check;
      if(PRINT_PATH_STEP_PROB)
	{
	  indiv_pois->push_back(temp_inf_check);
	}
      i++;
      //cout << ", loglike = " << probability_sum << endl;
    }
  return(probability_sum);
}

/*Actually now I think about it, this isn't quite right.  It should not be conditional poisson, but should
  count only the 0-0 transitions.  So if we have a single off state, we transition TO that state but we have 
  a count of 0 for the offlambda transitioner?*/

double CondPoissonLoglike(int x, double lambda)
{
  //Loglikelihood of poisson given greater than 0
  double val;
  val = log(gsl_ran_poisson_pdf(x, lambda));// - log(1 - gsl_ran_poisson_pdf(0,lambda));
  return(val);
}



double GetLogTransitionProbabilities(mach_dat_obj* data, curr_par_obj* current_values, global_par_obj* global_pars, int startval, int endval, vector<double>* indiv_trans)
{
  int i, end_trunc;
  double loglike;
  int not_infinity;
  double temp_inf_check;
  double rho, gamma, nu;
  int offcounter;
  
  if(endval >= current_values->Deathdate())
    {
      end_trunc = current_values->Deathdate();
      if(end_trunc == current_values->T)
	end_trunc--;
    }
  else
    end_trunc = endval;

  if(PRINT_PATH_STEP_PROB)
    {
      while(indiv_trans->size() > 0)
	indiv_trans->pop_back();
      indiv_trans->push_back(0.0);
    }

  loglike = 0;
  if(startval >= current_values->Birthdate())
    i = startval;
  else
    i = current_values->Birthdate();
  not_infinity = 1;
  offcounter = 0;

  while( (i <= (end_trunc-1)))
    {
      if(current_values->State(i) != 0)//we have a non-zero state
	{
	  if(offcounter > 0)
	    {
	      temp_inf_check = CondPoissonLoglike(offcounter-1, current_values->off_lambda);
	      loglike += temp_inf_check;
	      offcounter = 0;
	      if(PRINT_PATH_STEP_PROB)
		indiv_trans->operator[](indiv_trans->size() - 1) = temp_inf_check;
	    }
	  rho = current_values->rho[current_values->State(i)];
	  gamma = current_values->gamma[current_values->State(i)];
	  nu = current_values->nu[current_values->State(i)];
	  temp_inf_check = log(sinecurveprob_single(current_values->State(i), current_values->State(i + 1), i+1 - current_values->tstar, rho, gamma, nu, current_values));
      //cout << endl << "transitionprobs: cur_state = " << current_values->state[i] << ", next_state = " << current_values->state[i+1] << ", temp_check = "<< temp_inf_check; 
	  if(temp_inf_check == (-1.0)*INFINITY)
	    {
	      not_infinity = 0;
	      loglike = temp_inf_check;
	    }
	  if(not_infinity)
	    loglike += temp_inf_check;
	  if(PRINT_PATH_STEP_PROB)
	    {
	      indiv_trans->push_back(temp_inf_check);
	    }
	  if(!(current_values->Is_immune()))
            {
              if( (i+1) != current_values->Deathdate())
		temp_inf_check += log(global_pars->survival_rate[i]);
	      else
		temp_inf_check += log((1.0 - global_pars->survival_rate[i]));
	      if(temp_inf_check == (-1.0)*INFINITY)
		{
		  not_infinity = 0;
		  loglike = temp_inf_check;
		}
	      if(not_infinity)
		loglike += temp_inf_check;
	    }
	  // cout << ", loglike = " << loglike << endl;
	}
      else
	{
	  offcounter++;
	  if(PRINT_PATH_STEP_PROB)
	    indiv_trans->push_back(0.0);
	}
      i++;
    }
  if((offcounter > 0))
    {     
      temp_inf_check = CondPoissonLoglike(offcounter-1, current_values->off_lambda);
      if(not_infinity)
	loglike += temp_inf_check;
      if(PRINT_PATH_STEP_PROB)
	indiv_trans->operator[](indiv_trans->size() - 1) = temp_inf_check;
    }
  
  return(loglike);
  
}


double GetPathLogLikelihood(int startval, int endval, mach_dat_obj* data, curr_par_obj* current_values,global_par_obj* global_pars, int whichtype, vector<double>* indiv_trans, vector<double>* indiv_pois)
{
  //whichtype: 0 = full, 1 = path probabilities, 2 = poisson probabilities

  double loglike;
  double tempval_inf_check;
  loglike = 0;
  if( (whichtype == 0) || (whichtype == 2))
    {
      tempval_inf_check =  GetLogPoissonProbabilities(data, current_values, startval, endval, indiv_pois);
      if(tempval_inf_check != (-1.0*INFINITY))
	loglike += tempval_inf_check;
      else
	loglike = ( -1.0 * INFINITY);
    }
  if( ((whichtype == 0) || (whichtype == 1)) && (loglike != (-1.0*INFINITY))   ) 
    {
      tempval_inf_check = GetLogTransitionProbabilities(data, current_values, global_pars, startval, endval, indiv_trans);
      if(tempval_inf_check != (-1.0 * INFINITY))
	loglike += tempval_inf_check;
      else
	loglike = (-1.0 * INFINITY);
    }
  return(loglike);
}



double EtaPropH(int index, int startval, int endval, mach_dat_obj* data, curr_par_obj* current_values, global_par_obj* global_pars, vector<double>* indiv_trans, vector<double>* indiv_pois)
{
  double tpar;
  tpar = GetPathLogLikelihood(startval, endval, data, current_values, global_pars, 0, indiv_trans, indiv_pois);
  return(tpar);
}


double LogAplusB(double loga, double logb)
{
  double minval;
  double maxval;
  double ret;
  
  if((loga*logb != 0) && (loga != (-1.0*INFINITY)) && (logb != (-1.0*INFINITY)))
    {
      if(loga > logb)
	{
	  maxval = loga;
	  minval = logb;
	}
      else
	{
	  maxval = logb;
	  minval = loga;
	}
      ret = maxval + gsl_log1p(exp(minval - maxval));
      //cout << "loga logb = "<< loga << " " << logb << ": non-zero, non-inf, exp(min - max) = " << exp(minval - maxval) << endl;
    }
  else
    {
      if((loga == 0) || (loga == (-1.0*INFINITY)))
	ret = logb;
      else
	ret = loga;
    }
  return(ret);
}

double BetaDistAlpha(double mean, double temperature)
{
  double val;
  if(  (mean < 0) || (mean > 1) || (temperature < 1)  )
    {
      cout << "Error: BetaDistAlpha(double mean, double temperature): Invalid mean = " << mean << " temperature = " << temperature << " specification: constraints are {0 <= mean <= 1,  1 <= temperature}" <<  endl;
      val = -1.0;
    }
  else
    {
      val = mean*(temperature - 1.0);
    }
  return val;
}

double BetaDistBeta(double mean, double temperature)
{
  double val;
  if(  (mean < 0) || (mean > 1) || (temperature < 1)  )
    {
      cout << "Error: BetaDistAlpha(double mean, double temperature): Invalid mean = " << mean << " temperature = " << temperature << " specification: constraints are {0 <= mean <= 1,  1 <= temperature}" <<  endl;
      val = -1.0;
    }
  else
    {
      val = temperature - (mean*(temperature -1)) - 1.0;
    }
  return val;
}


double ScaledGammaLogConstant(double lower, double upper, double alpha, double beta)
{
  double val;
  val = gsl_cdf_gamma_P(upper, alpha, beta) - gsl_cdf_gamma_P(lower, alpha, beta);
  if(val < GAMMA_SCALE_MIN)
    val = GAMMA_SCALE_MIN;
  return(log(val));   
}

double LogSpikeProb(double spikemult, short int state1, short int state2)
{
  double val;
  //q(state1 | state2) under spike multiplier spikemult
  if(state1 == 1)
    {
      val = 1.0-spikemult;
    }
  else
    {
      if(state2 == 1)
	val = 0.5;
      else
	val = spikemult;
    }
  return log(val);
}

void StateRunChanges(curr_par_obj* current_values, curr_par_obj* proposed_values, int index, int beginval, int endval, int* add)
{
  short int oldstate, newstate;
  short int pristate0;
  short int poststate0;
  int addsub = 1;
  pristate0 = 0;
  poststate0 = 0;
  add[0] = 0;
  add[1] = 0;
  oldstate = current_values->State(index);
  newstate = proposed_values->State(index);
  if( ((oldstate == 0) || (newstate == 0)) && (oldstate != newstate) )
    {
      if(beginval < (index-1))
	pristate0 = (current_values->State(index-1) == 0);
      if(endval > (index+1))
	poststate0 = (current_values->State(index+1) == 0);
      if(oldstate == 0)
	addsub = -1;

      add[0] = (addsub)*(pristate0 + poststate0);
      add[1] = (-1*addsub)*(pristate0 + poststate0 -1);
    }
  return;
}

double LogGamma(double x)
{
  return gsl_sf_lngamma(x);
}

double LogPriorRatioScaledGamma(double s1, double s2, double s, double minval, double maxval, double pri_k, double pri_theta)
{
  //default is merge, multiply by -1 to get split
  double ratio;
  double loglike;
  double normconst = 1.0;
  double upper_cdf, lower_cdf;
  ratio = s/(s1*s2);
  //avoid -infinities

  lower_cdf = gsl_cdf_gamma_P(minval, pri_k, pri_theta);
  if(maxval < 0)
    {
      //infinite
      upper_cdf = 1.0;
    }
  else
    {
      upper_cdf = gsl_cdf_gamma_P(maxval, pri_k, pri_theta);
    }
  normconst = upper_cdf - lower_cdf;
  if(normconst < GAMMA_SCALE_MIN)
    {
      //it really shouldn't be
      normconst = GAMMA_SCALE_MIN;
    }
  if(pri_theta < MIN_LOG_P)
    pri_theta = MIN_LOG_P;
  if(pri_k < MIN_LOG_P)
    pri_k = MIN_LOG_P;
  if(ratio < MIN_LOG_P)
    ratio = MIN_LOG_P;  
  loglike = log(normconst) + pri_k*log(pri_theta) + LogGamma(pri_k) + (pri_k - 1.0)*log(ratio) + (s1 + s2 - s)/pri_theta;
  return loglike;
}

double LogPriorRatioScaledBeta(double s1, double s2, double s, double minval, double maxval, double pri_alpha, double pri_beta)
{
  double ratio_alpha, ratio_beta;
  double loglike;
  ratio_alpha = ((maxval - minval)*(s - minval))/( (s1-minval)*(s2-minval));
  if(ratio_alpha < MIN_LOG_P)
    ratio_alpha = MIN_LOG_P;
  if(ratio_alpha > MAX_LOG_P)
    ratio_alpha = MAX_LOG_P;
  
  ratio_beta = ((maxval - minval)*(maxval - s))/( (maxval - s1)*(maxval - s2));
  if(ratio_beta < MIN_LOG_P)
    ratio_beta = MIN_LOG_P;
  if(ratio_beta > MAX_LOG_P)
    ratio_beta = MAX_LOG_P;

  loglike = log(maxval - minval) + LogGamma(pri_alpha) + LogGamma(pri_beta) - LogGamma(pri_alpha + pri_beta) + (pri_alpha - 1.0)*log(ratio_alpha) + (pri_beta - 1.0)*log(ratio_beta);
  return loglike;
}

int FirstLastRun(curr_par_obj* current_values, int which, int use_absolute)
{
  //which = 1 for first run from birth, -1 for last run from death
  //if use_absolute == True, get first nonoff from 0, last nonoff from T
  int found;
  int index;
  int is_at_bound = 0;
  found = 0;
  if(which == 1)
    {
      if(use_absolute)
	index = 0;
      else
	index = current_values->Birthdate();
    }
  else
    {
      if(use_absolute)
	{
	  index = current_values->T - 1;
	}
      else
	{
	  index = current_values->Deathdate();
	  if(index >= current_values->T)
	    index = current_values->T -1;
	}
    }
  if(current_values->State(index) != 0)
    found = 1;
  else 
    {
      while((!found) && (!(is_at_bound))) 
	{
	  index += which;
	  if(((which==1)&&(index==current_values->Deathdate()))||((which==-1)&&(index==current_values->Birthdate()))) 
	    {
	      is_at_bound = 1;
	    }
	  if(current_values->State(index) != 0)
	    {
	      found = 1;
	    }
	}
    }
  if( (is_at_bound) && !found)
    {
      cout << "Error in FirstLastRun: Bound reached but first (last) nonoff state not found" << endl;
    }
  return index;
}


int GetBDIIndex(short int birth, short int death, short int immune, int max_death)
{
  int subind = 0;
  int index;
  index = 4*(birth != 0) + 2*(death == max_death) + (immune == 1);
  if(index > 0)
    subind++;
  if(index > 5)
    subind++;
  index = index - subind;
  return(index);
}

void BDIFromIndex(int index, int last_nonoff, int max_death, short int* birth, short int* death, short int* immune)
{
  if(index < 3)
    {
      (*birth) = 0;
    }
  else
    {
      (*birth) = 284;
    }
  if( (index == 0) || (index == 3))
    {
      (*death) = last_nonoff;
    }
  else
    {
      (*death) = max_death;
    }
  if( (index == 2) || (index == 5))
    {
      (*immune) = 1;
    }
  else
    {
      (*immune) = 0;
    }
  return;
}

double LogMoveRatio(vector<double>& moveprobs, int proposed_type, int current_type)
{
  //move from 
  double val = 0;
  if( (proposed_type > moveprobs.size()) || (current_type > moveprobs.size()))
    cout << "Error in LogMoveRatio: Types not found" << endl;
  else
    val = log(moveprobs[current_type]) - log(moveprobs[proposed_type]);
  return val;
}

int SimpleSample(vector<double>& prob, DistributionList* dists, double probsum)
{
  int index = 0;
  double sumval = 0;
  double cursum;
  int i;
  double val;
  if(probsum > 0)
    sumval = probsum;
  else
    {
      for(i = 0; i < prob.size(); i++)
	{
	  sumval += prob[i];
	}
    }

  dists->Simulate(3, 1, &val, 0.0, sumval);
  cursum = prob[0];
  while((cursum < val) && (index+1 < prob.size()))
    {
      index++;
      cursum += prob[index];
    }
  return index;
}


void GetPlusMinusState(int state, curr_par_obj* current_values, double& plusstatecount, double& minusstatecount)
{
  int i, end_trunc;
  double pvals[3];
  int plusgammastate;
  int minusgammastate;

  plusstatecount = 0.0;
  minusstatecount = 0.0;
  end_trunc = current_values->Deathdate();
  if(end_trunc >= current_values->T)
    {
      end_trunc = current_values->T -1;
    }
  pvals[0] = current_values->psi_parameterization[state][0][1];
  pvals[1] = current_values->psi_parameterization[state][1][1];
  pvals[2] = current_values->psi_parameterization[state][2][1];
  
  for(i = 0; i < 3; i++)
    {
      if(pvals[i] == 1)
	plusgammastate = i;
      if(pvals[i] == -1)
	minusgammastate = i;
    }
  
  for(i = current_values->Birthdate() + 1; i <= end_trunc; i++)
    {
      if(current_values->State(i-1) == state)
	{
	  if(current_values->State(i) == plusgammastate)
	    plusstatecount+= 1.0;
	  else
	    if(current_values->State(i) == minusgammastate)
	      minusstatecount+= 1.0;
	  
	}
      
    }
  
}


void RhoNuSuccessesAndTotalCounts(short int state, curr_par_obj* current_values, int mod_index, int& successes, int& totals)
{
  //assuming successes and total counts are both of length 24
  int i;
  int bdate, ddate;
  int cur_state, next_state;
  int psi_val;
 
  successes = 0;
  totals = 0;
  bdate = current_values->Birthdate();
  ddate = current_values->Deathdate();
  cur_state = current_values->State(bdate);
  if(ddate >= current_values->T)
    ddate = current_values->T - 1;
  
  for(i = bdate; i < ddate-1; i++)
    {
      next_state = current_values->State(i+1);
      if(cur_state == state)
	{
	  //if we are at the relavent state
	  if( ((i - current_values->tstar)%24) == mod_index )
	    {
	      psi_val = current_values->psi_parameterization[cur_state][next_state][0];
	      if(psi_val != 0)
		{ 
		  totals += 1;
		}
	      if(psi_val == 1)
		{
		  successes += 1;
		}	    
	    }
	}
      cur_state = next_state;      
    }
  
}


void ActiveOverlap(curr_par_obj* par1, curr_par_obj* par2, double& A1, double& A2, double& A1U2, double& A1I2)
{
  A1 = 0;
  A2 = 0;
  A1U2 = 0;
  A1I2 = 0;
  int s1, s2;
  int i;
  for(i = 0; i < par1->T; i++)
    {
      s1 = par1->State(i);
      s2 = par2->State(i);
      if(s1 > 0)
	{
	  A1++;
	}
      if(s2 > 0)
	{
	  A2++;
	}
      if( (s1>0) || (s2> 0))
	{
	  A1U2++;
	}
      if( (s1>0) && (s2> 0))
	{
	  A1I2++;
	}
	  
    }
  
}

double ScaledBetaPriorMode(double l, double u, double alpha, double beta)
{
  double m;
  if( (alpha < 1) && (beta < 1))
    {
      m = l - 1;
    }
  else
    {
      if(alpha < 1)
	m = l + BETA_MIN;
      else
	{
	  if(beta < 1)
	    {
	      m = u - BETA_MIN;
	    }
	  else
	    {
	      m = ((alpha - 1)*( u - l)/(alpha + beta - 2)) + l;
	    }
	}
    }
  
  return m;
}

void ScaledBetaAlphaBeta(double l, double u, double mode, double samplesize, double& alpha, double &beta)
{
  double tempalpha, tempbeta;
  double tempmode;
  if(mode < l)
    mode = l + BETA_MIN;
  if(mode > u) 
    mode = u - BETA_MIN;

  tempalpha = ((mode + l)*(samplesize-2)/(u - l)) + 1.0;
  tempbeta = samplesize - 1.0 - (mode*(samplesize-2.0)/(u - l));

  

  alpha = tempalpha;
  beta = tempbeta;
 
}

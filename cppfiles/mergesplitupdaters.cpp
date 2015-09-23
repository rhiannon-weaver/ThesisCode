#include "../mergesplitupdaters.h"



/*===================Helper Functions for Common Merge/Split Actions====================*/

/*============Asymmetric Scaled Average (2D-2D)=====================*/
AsymmetricAverage::AsymmetricAverage(double lb, double ub, double prop_alpha, double prop_beta, int is_u_infinity, int is_l_infinity, double assignment_prob_unif_weight, double assignment_prob_power, DistributionList* distributions )
{
  l = lb;
  u = ub;
  l_inf = is_l_infinity;
  u_inf = is_u_infinity;
  dists = distributions;
  assign_alpha = assignment_prob_unif_weight;
  assign_k = assignment_prob_power;
  u1_alpha = prop_alpha;
  u1_beta = prop_beta;
}

void AsymmetricAverage::MergePropose(double q1, double q2, double &q)
{
  double qmin, qmax,qtemp, u1;
  if(q1 < q2)
    {
      qmin = q1;
      qmax = q2;
    }
  else
    {
      qmin = q2;
      qmax = q1;
    }
  ginv(qmin, qmax, qtemp, u1);
  q = qtemp;
}

void AsymmetricAverage::SplitPropose(double q, double &q1, double &q2, double assignment_prob)
{
  double qmin, qmax;
  int assign_alpha;
  double u1;  
  dists->Simulate(7, 1, &u1, 0.0, u-q, u1_alpha, u1_beta);
  gfun(q, u1, qmin, qmax);
  if(qmin != qmax)
    {
      dists->Simulate(3,1,&assign_alpha, 0.0, 1.0);
      if(log(assign_alpha) < AdjustAssignment(assignment_prob,0))
	{
	  q1 = qmax;
	  q2 = qmin;
	}
      else
	{
	  q2 = qmax;
	  q1 = qmin;
	}
    }
  else
    {
      q1 = qmax;
      q2 = qmin;
    }
}

double AsymmetricAverage::MergeLogRatio(double q1, double q2, double assignment_prob)
{
  double q, u1, det;
  double logratio;
  double logalpha = 0.0;
  if(q1 < q2)
    {
      logalpha = AdjustAssignment(assignment_prob, 1);
    }
  else
    {
      if(q1 > q2)
	logalpha = AdjustAssignment(assignment_prob, 0);
    }
  ginv(q1, q2, q, u1);
  det = JacobianDqu1Dq1q2(q1, q2);
 
  logratio = dists->LogDensity(7,u1, l, u, u1_alpha, u1_beta) + logalpha + log(det);
  return(logratio);
}

double AsymmetricAverage::SplitLogRatio(double q1, double q2, double assignment_prob)
{
  double q, u1, det;
  double u_alpha, u_beta;
  double logratio;
  double logalpha = 0.0;
  if(q1 < q2)
    {
      logalpha = AdjustAssignment(assignment_prob, 1);
    }
  else
    {
      if(q1 > q2)
	logalpha = AdjustAssignment(assignment_prob, 0);
    }
  ginv(q1, q2, q, u1);
  det = JacobianDq1q2Dqu1(q, u1);
 
  logratio = log(det) - dists->LogDensity(7, u1, l, u, u1_alpha, u1_beta);
  return(logratio);
}
 
void AsymmetricAverage::gfun(double q, double u1, double& qmin, double& qmax)
{
  qmax = q + u1;
  qmin = q - (u1/(u-q))*(q - l);
}
 
void AsymmetricAverage::ginv(double q1, double q2, double& q, double& u1)
{
  double qmin, qmax;
  if(q1 < q2)
    {
      qmin = q1;
      qmax = q2;
    }
  else
    {
      qmin = q2;
      qmax = q1;
    }
  q = ( (qmin * u) - (qmax*l))/( (u-qmax) + (qmin-l) );
  u1 = qmax - ( (qmin*u) - (qmax*l))/( (u-qmax) + (qmin-l) );
}
 
double AsymmetricAverage::JacobianDqu1Dq1q2(double q1, double q2)
{
  double qmin, qmax;
  double det;
  double denom;
  double dq_dqmin;
  if(q1 < q2)
    {
      qmin = q1;
      qmax = q2;
    }
  else
    {
      qmin = q2;
      qmax = q1;
    }

  denom = (u - qmax) + (qmin - l);
  dq_dqmin = (qmax*u - qmin*l)/(denom*denom) + (l)/denom;
  det = dq_dqmin;
  if(det < 0)
    det = -1.0*det;
  return(det);
}

double AsymmetricAverage::JacobianDq1q2Dqu1(double q, double u1)
{
  double a, b, c, d, det;
  a = 1.0;
  b = 1.0;
  c = 1.0 - (u1/((u-q)*(q-l))) +  ((u1*(q-l))/((u-q)*(u-q)*(u-q)));
  d = (q-l)/(u-q);
  det = a*d - b*c;
  if(det < 0)
    det = -1.0*det;
  return(det);
}

double AsymmetricAverage::AdjustAssignment(double assignment_prob, double invert)
{
  // a = ((alpha)*p + (1-alpha)*0.5)^k
  double a;
  a = (assign_alpha)*(assignment_prob) + (1.0 - assign_alpha)*(0.5);
  if(!invert)
    a = assign_k * log(a);
  else
    a = log(1.0 - exp(assign_k*log(a)));
  return a;
}

/*===========End Asymmetric Average============*/


/*=======Bounded Rescaled Sampler (3D-3D)=====*/

BoundedRescaledSampler::BoundedRescaledSampler(double lb, double ub, double alphaplusbeta, DistributionList* dists)
{
  l = lb;
  u = ub;
  samplesize = alphaplusbeta;
  distributions = dists;
}

void BoundedRescaledSampler::MergePropose(double q1, double q2, double & q,  vector<double>& proposal_pars,  vector<Machine*>& current_machines, vector<Machine*>& proposed_machines)
{
  double meanval, alpha, beta, res;
  gfun2(q1, q2, meanval,proposal_pars, current_machines, proposed_machines);
  GetAlphaBeta(meanval, alpha, beta);
  //ScaledBetaAlphaBeta(l, u, meanval, samplesize, alpha, beta);
  distributions->Simulate(7, 1, &res, l, u, alpha, beta);
  if(res >= u)
    res = u - BETA_MIN;
  if(res <= l)
    res = l + BETA_MIN;
  q = res;
}

void BoundedRescaledSampler::SplitPropose(double q, double & q1, double & q2, vector<double>& proposal_pars,  vector<Machine*>& current_machines, vector<Machine*>& proposed_machines)
{
  double meanval1, meanval2, alpha, beta, res1, res2;
  gfun1(q, meanval1, meanval2, proposal_pars, current_machines, proposed_machines);
  //cout << "mean transform: q = "<< q << " -> " << meanval1 << "," << meanval2 << endl;
  GetAlphaBeta(meanval1, alpha, beta);
  //ScaledBetaAlphaBeta(l, u, meanval1, samplesize, alpha, beta);
  distributions->Simulate(7, 1, &res1, l, u, alpha, beta);
  if(res1 >= u)
    res1 = u - BETA_MIN;
  if(res1 <= l)
    res1 = l + BETA_MIN;
  GetAlphaBeta(meanval2, alpha, beta);
  //ScaledBetaAlphaBeta(l, u, meanval2, samplesize, alpha, beta);
  distributions->Simulate(7, 1, &res2, l, u, alpha, beta);
  if(res2 >= u)
    res2 = u - BETA_MIN;
  if(res2 <= l)
    res2 = l + BETA_MIN;

  q1 = res1;
  q2 = res2;
}

void BoundedRescaledSampler::GetAlphaBeta(double q, double &alpha, double& beta)
{
  double aval, bval;
  double meanval_01 = (q - l)/(u - l);
  alpha = meanval_01*samplesize;
  beta = (1-meanval_01)*samplesize;
  //cout << "beta transformation: mu = "<< meanval_01 << " ss = "<<samplesize << " -> alpha = " << alpha << " beta = " << beta << endl;
}

double BoundedRescaledSampler::LogRatio(double q1, double q2, double q12, vector<double>& proposal_pars,   vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, int mergesplit)
{
  //default merge:  d(1) + d(2) - d(12) 
  double meanval_12,meanval_1, meanval_2;
  double alpha_12, beta_12, alpha_1, beta_1, alpha_2, beta_2;
  double loglike;
  gfun1(q12,meanval_1, meanval_2, proposal_pars, current_machines, proposed_machines);
  gfun2(q1, q2, meanval_12,proposal_pars, current_machines, proposed_machines);
  GetAlphaBeta(meanval_12, alpha_12, beta_12);
  //ScaledBetaAlphaBeta(l, u, meanval_12, samplesize, alpha_12, beta_12);
  GetAlphaBeta(meanval_1, alpha_1, beta_1);
  //ScaledBetaAlphaBeta(l, u, meanval_1, samplesize, alpha_1, beta_1);
  GetAlphaBeta(meanval_2, alpha_2, beta_2);
  //ScaledBetaAlphaBeta(l, u, meanval_2, samplesize, alpha_2, beta_2);
 
  loglike = distributions->LogDensity(7, q1, l, u, alpha_1, beta_1) + distributions->LogDensity(7, q2, l, u, alpha_1, beta_1) - distributions->LogDensity(7, q12, l, u, alpha_12, beta_12);
  
  return( ((double)(mergesplit))*loglike);

}

void BoundedRescaledSampler::gfun1(double q, double & q1, double & q2,  vector<double>& proposal_pars,  vector<Machine*>& current_machines, vector<Machine*>& proposed_machines)
{
  q1 = q;
  q2 = q;

  if(q1 < l)
	  q1 = l + BETA_MIN;
  if(q2 < l)
	  q2 = l + BETA_MIN;
	
  if(q1 > u)
	  q1 = u - BETA_MIN;
  if(q2 > u)
	  q2 = u - BETA_MIN;
	
}

void BoundedRescaledSampler::gfun2(double q1, double q2, double & q,  vector<double>& proposal_pars,  vector<Machine*>& current_machines, vector<Machine*>& proposed_machines)
{
  q = 0.5*(q1 + q2);
  if(q < l)
	  q = l + BETA_MIN;
  if(q > u)
	  q = u - BETA_MIN;
}

double BoundedRescaledSampler::MergeLogRatio(double q1, double q2, double q12,  vector<double>& proposal_pars,  vector<Machine*>& current_machines, vector<Machine*>& proposed_machines)
{
  return(LogRatio( q1, q2, q12,proposal_pars, current_machines, proposed_machines, 1));
}

double BoundedRescaledSampler::SplitLogRatio(double q1, double q2, double q12,  vector<double>& proposal_pars,  vector<Machine*>& current_machines, vector<Machine*>& proposed_machines)
{
  return(LogRatio( q1, q2, q12, proposal_pars, current_machines, proposed_machines, -1));
}





/*================End Helper Functions for Common Merge/Split Actions===================*/





/*==============================Base Class: MergeSplitUpdater============================*/

MergeSplitUpdater::MergeSplitUpdater(vector<double>&hpar, DistributionList* dists, ParameterUpdater* parameter_updater)
{
  int i;
  for(i = 0; i < hpar.size(); i++)
    proposal_par.push_back(hpar[i]);
  for(i = 0; i < parameter_updater->num_hyperpars; i++)
    prior_hyperpar.push_back(parameter_updater->hyperpar[i]);
  par_name = "MergeSplit_" + parameter_updater->par_name;
  distributions = dists;
}

MergeSplitUpdater::MergeSplitUpdater(vector<double>&hpar, DistributionList* dists)
{
  int i;
  for(i = 0; i < hpar.size(); i++)
    proposal_par.push_back(hpar[i]);
  distributions = dists;
}

double MergeSplitUpdater::MergeSplitStep(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars, int print_to_screen)
{
  double loglike;
  if(current_machines.size() == 1)
    {
      SplitPropose(networks, current_machines, proposed_machines, global_pars);
      loglike = SplitLogRatio(networks, current_machines, proposed_machines, global_pars);
    }
  else
    {
      if(current_machines.size() == 2)
	{
	  MergePropose(networks, current_machines, proposed_machines, global_pars);
	  loglike = MergeLogRatio(networks, current_machines, proposed_machines, global_pars);
	}
      else
	{
	  cout << "Error, current machines size = " << current_machines.size();
	}
    }
  if(print_to_screen)
    PrintToScreen(networks, current_machines, proposed_machines, loglike);
  return(loglike);
}
 
void MergeSplitUpdater::PrintToScreen(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, double loglike)
{
  cout << par_name <<": ";
  if(current_machines.size() == 1)
    {
      cout << "Split Step";
    }
  else
    {
      if(current_machines.size() == 2)
	cout << "Merge Step";
      else
	cout << "Error, current machines size = " << current_machines.size();
    }
  cout << "| current= "; 
  PrintValue(current_machines, cout);
  cout << " | proposed= ";
  PrintValue(proposed_machines, cout);
  cout << " | log-likelihood=" << loglike << endl;
}

void MergeSplitUpdater::PrintHyperpars(ostream& str)
{
  int i;
  str << par_name << " (";
  for(i = 0; i < prior_hyperpar.size(); i++)
    str << setprecision(8)<<prior_hyperpar[i] << " " ;
  str << ") (";
  for(i = 0; i < proposal_par.size(); i++)
    str <<  setprecision(8)<<proposal_par[i] << " ";
  str << ")" << endl;
}

/*===========================End Base Class: MergeSplitUpdater============================*/


/*==========================Derived Classes: Individual Parameters========================*/


/*=====================================State=================================*/

StateMergeSplitUpdater::StateMergeSplitUpdater(vector<double>& hpar1, DistributionList* dists, ParameterUpdater* par):MergeSplitUpdater(hpar1, dists, par)
 {
   y_active_stats = new YTotalActiveStats();

 }

void StateMergeSplitUpdater::MergePropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
{
  int i;
  short int s, s1, s2;
  for(i = 0; i< 3;  i++)
    {
      s_count[i] = 0;
      s1_count[i] = 0;
      s2_count[i] = 0;
    }
  y_active_stats->Initialize(current_machines);
  

  for(i = 0; i < current_machines[0]->parameters->T; i++)
    {
      s1 = current_machines[0]->parameters->State(i);
      s2 = current_machines[1]->parameters->State(i);
      s1_count[s1]++;
      s2_count[s2]++;
      s = MergeState(s1, s2);
      s_count[s]++;
      proposed_machines[0]->parameters->SetState(i, s);
    }
  
}

void StateMergeSplitUpdater::SplitPropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
{
  int i;
  double pri_lagstate1, pri_lagstate2;
  double pri_offrun1, pri_offrun2;
  double prob_both_on;
  short int s, s1, s2;
  for(i = 0; i< 3;  i++)
    {
      s_count[i] = 0;
      s1_count[i] = 0;
      s2_count[i] = 0;
    }
  y_active_stats->Initialize(current_machines);
  log_split_prob = 0.0;
  SetLagstateAndOffrun(proposed_machines, -1, pri_lagstate1, pri_lagstate2, pri_offrun1, pri_offrun2);
  for(i = 0; i < current_machines[0]->parameters->T; i++)
    {
      s = current_machines[0]->parameters->State(i);
      s_count[s]++;							
      prob_both_on = GetProbabilityBothOn(current_machines, i);

      log_split_prob += SplitState(s, prob_both_on,pri_lagstate1, pri_lagstate2, pri_offrun1, pri_offrun2, s1, s2);
      
      s1_count[s1]++;
      s2_count[s2]++;
      
      proposed_machines[0]->parameters->SetState(i, s1);
      proposed_machines[1]->parameters->SetState(i, s2);
      SetLagstateAndOffrun(proposed_machines, i, pri_lagstate1, pri_lagstate2, pri_offrun1, pri_offrun2);
    }
}

double StateMergeSplitUpdater::GetProbabilityBothOn(vector<Machine*>& machines, int i)
{
  int y, z;
  y = machines[0]->data->GetTotalCount(i);
  for(z = 1; z < machines.size(); z++)
    y += machines[z]->data->GetTotalCount(i);
  return( y_active_stats->YCDF(y)  );
}

void StateMergeSplitUpdater::SetLagstateAndOffrun(vector<Machine*>& machines, int i, double& pri_lagstate1, double& pri_lagstate2, double& pri_offrun1, double& pri_offrun2)
{
  
  if(i == -1)
    {
        pri_lagstate1 = LAG_MAX;
	pri_lagstate2 = LAG_MAX;
	pri_offrun1 = 0.0;
	pri_offrun2 = 0.0;
    }
  else
    {
      pri_lagstate1 = machines[0]->parameters->Lagstate(i);
      if(pri_lagstate1 > LAG_MAX)
	pri_lagstate1 = LAG_MAX;
      if(pri_lagstate1 == -1)
	{
	  pri_offrun1 += 1.0;
	  if(pri_offrun1 > LAG_MAX)
	    pri_offrun1 = LAG_MAX;
	}
      else
	pri_offrun1 = 0.0;
      pri_lagstate2 = machines[1]->parameters->Lagstate(i);
      if(pri_lagstate2 > LAG_MAX)
	pri_lagstate2 = LAG_MAX;
      if(pri_lagstate1 == -1)
	{
	  pri_offrun2 += 1.0;
	  if(pri_offrun2 > LAG_MAX)
	    pri_offrun2 = LAG_MAX;
	}
      else
	pri_offrun2 = 0.0;
    }
}

double StateMergeSplitUpdater::GetProb1Off(double pri_lagstate1, double pri_lagstate2, double pri_offrun1, double pri_offrun2)
{
  double exponent = 1.0;
  int reverse = 0;
  double prob;

  if( (pri_lagstate1 >= 0) && (pri_lagstate2 >=0))
    {
      exponent = (pri_lagstate1 + 1.0)/(pri_lagstate2 + 1.0);   
    }
  else
    {
      if((pri_lagstate1 == -1) && (pri_lagstate2 == -1))
	{
	  exponent = (pri_offrun1 + 1.0)/(pri_offrun2 + 1.0);
	}
      else
	{
	  if(pri_lagstate1 == -1)
	    {
	      exponent = pri_offrun1 + 1;
	      reverse  = 1;
	    }
	  else
	    {
	      exponent = pri_offrun2 + 1;
	      reverse = 1;
	    }
	}
    }
  
  prob = exp(exponent*log(0.5));
  if(reverse)
    prob = 1.0 - prob;
  return prob;
}


double StateMergeSplitUpdater::MergeLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
{
  //prod_t g(1->2)/g(2->1): ie sum_t log(g(1->2)) - log(g(2->1))
  double pri_lagstate1, pri_lagstate2, pri_offrun1, pri_offrun2;
  short int s1, s2, s;
  double prob_both_on;
  int i;
  double loglike = 0.0;
  
  SetLagstateAndOffrun(current_machines, -1, pri_lagstate1, pri_lagstate2, pri_offrun1, pri_offrun2);
  for(i = 0; i < current_machines[0]->parameters->T; i++)
    {
      prob_both_on = GetProbabilityBothOn(current_machines, i);
      s1 = current_machines[0]->parameters->State(i); 
      s2 = current_machines[1]->parameters->State(i);
      s = proposed_machines[0]->parameters->State(i);
      
      loglike += log(SplitStateProb(s, s1, s2, prob_both_on, pri_lagstate1, pri_lagstate2, pri_offrun1, pri_offrun2)) - log(MergeStateProb(s1, s2, s)); 
      SetLagstateAndOffrun(current_machines, i, pri_lagstate1, pri_lagstate2, pri_offrun1, pri_offrun2);
    }
  return loglike;
}


double StateMergeSplitUpdater::SplitLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
{
  //prod_t g(2->1)/g(1->2): ie sum_t log(g(2->1)) - log(g(1->2))
  short int s1, s2, s;
  int i;
  double loglike = 0.0 - log_split_prob;
  
  for(i = 0; i < current_machines[0]->parameters->T; i++)
    {
      s1 = proposed_machines[0]->parameters->State(i);
      s2 = proposed_machines[1]->parameters->State(i);
      s = current_machines[0]->parameters->State(i);
      loglike += log(MergeStateProb(s1, s2, s));
    }
  return loglike;
}

void StateMergeSplitUpdater::PrintValue(vector<Machine*>& machines, ostream& str)
{
  if(machines.size() == 1)
    {
      str << machines[0]->uid << ": (" << s_count[0] << " off, " << s_count[1] << " spike, " << s_count[2] << " decay) ";
    }
  else
    {
      str << machines[0]->uid << ": (" << s1_count[0] << " off, " << s1_count[1] << " spike, " << s1_count[2] << " decay) ";
      str << machines[1]->uid << ": (" << s2_count[0] << " off, " << s2_count[1] << " spike, " << s2_count[2] << " decay) ";
    }
}

short int StateMergeSplitUpdater::MergeState(short int s1, short int s2)
{
  //0 1 2 : off spike decay
  short int s;
  if( (s1+s2) == 3)
    s = 1;
  else
    {
      if(s2 > s1)
	s = s2;
      else
	s = s1;
    }
  return s;
}

double StateMergeSplitUpdater::SplitState(short int s, double p_both_on, double pri_lagstate1, double pri_lagstate2, double pri_offrun1, double pri_offrun2, short int& s1, short int& s2)
{
  //0 1 2 : off spike decay
  vector<short int> possibilities;
  vector<double> probabilities;
  possibilities.clear();
  probabilities.clear();
  int val;
  double prob;
  int choice;

  switch(s){
  case 0:
    possibilities.push_back(0);   //0 0
    probabilities.push_back(1.0);
    break;
  case 2:
    possibilities.push_back(2);   //0 2
    probabilities.push_back(SplitStateProb(2,0,2,p_both_on, pri_lagstate1, pri_lagstate2, pri_offrun1, pri_offrun2 ));
    possibilities.push_back(6);   //2 0
    probabilities.push_back(SplitStateProb(2,2,0,p_both_on, pri_lagstate1, pri_lagstate2, pri_offrun1, pri_offrun2 ));
    possibilities.push_back(8);   //2 2
    probabilities.push_back(SplitStateProb(2,2,2,p_both_on, pri_lagstate1, pri_lagstate2, pri_offrun1, pri_offrun2 ));
    break;
  case 1:
    possibilities.push_back(1);  //0 1
    probabilities.push_back(SplitStateProb(1,0,1,p_both_on, pri_lagstate1, pri_lagstate2, pri_offrun1, pri_offrun2 ));
    possibilities.push_back(7);  //2 1
    probabilities.push_back(SplitStateProb(1,2,1,p_both_on, pri_lagstate1, pri_lagstate2, pri_offrun1, pri_offrun2 ));
    possibilities.push_back(3);  //1 0
    probabilities.push_back(SplitStateProb(1,1,0,p_both_on, pri_lagstate1, pri_lagstate2, pri_offrun1, pri_offrun2 ));
    possibilities.push_back(5);  //1 2
    probabilities.push_back(SplitStateProb(1,1,2,p_both_on, pri_lagstate1, pri_lagstate2, pri_offrun1, pri_offrun2 ));
    possibilities.push_back(4);  //1 1
    probabilities.push_back(SplitStateProb(1,1,1,p_both_on, pri_lagstate1, pri_lagstate2, pri_offrun1, pri_offrun2 ));
    break;
  }

  val = SimpleSample(probabilities, distributions, 1.0);
  //  distributions->Simulate(3, 1, &val, 0.0, (double)(possibilities.size()));
  choice = possibilities[val];
  prob = probabilities[val];
  s1 = choice/3;
  s2 = choice%3;
  possibilities.clear();
  probabilities.clear();
  return prob;
}

double StateMergeSplitUpdater::MergeStateProb(short int s1, short int s2, short int s)
{
  double val = 1.0;
  if(s == 0)
    {
      if( !( (s1==0) && (s2==0) ) )
	val = 0.0;
    }
  if(s == 1)
    {
      if(!((s1 ==1) || (s2 ==1) ))
	val = 0.0;
    }
  if(s == 2)
    {
      if(s1 != 2)
	{
	  if(s2 != 2)
	    val = 0.0;
	}
      else
	{
	  if(s2 == 1)
	    val = 0.0;
	}
    }
  if(val == 0.0)
    {
      cout << "Impossible merge probability: "<< s1<<s2 << "->" << s << endl;
    }
  return val;
}



double StateMergeSplitUpdater::SplitStateProb(short int s, short int s1, short int s2, double prob_both_on, double pri_lagstate1, double pri_lagstate2, double pri_offrun1, double pri_offrun2)
{
  double val = 1.0;
  double prob_one_off = GetProb1Off(pri_lagstate1, pri_lagstate2, pri_offrun1, pri_offrun2);
  if(s == 0)
    {
      if( (s1 == 0) && (s2 == 0))
	val = 1.0;
      else
	val = 0.0;
    }
  if(s == 1)
    {
      if( (s1 == 1) && (s2 == 1))
	{
	  val *= prob_both_on;
	}
      else
	{
	  val *= (1.0 - prob_both_on);
	  if( ((s1 == 1) && (s2 == 2)) || ( (s1 == 2) && (s2 == 1)))
	    {
	      val *= prob_both_on;
	      if(s1 == 1)
		val *= 1.0-prob_one_off;
	      else
		val *= prob_one_off;
	    }
	  else
	    {
	      if( ((s1 == 1) && (s2 == 0)) || ( (s1 == 0) && (s2 ==1)))
		{
		  val *= (1.0 - prob_both_on);
		  if(s1 == 1)
		    val *= 1.0 - prob_one_off;
		  else
		    val *= prob_one_off;
		}
	      else
		val = 0.0;
	    }
	  
	}
    }
  if(s == 2)
    {
      if( (s1 == 2) && (s2 == 2))
	{
	  val *= prob_both_on;
	}
      else
	{
	  if(  ((s1 == 2)&&(s2 == 0)) || ((s1 == 0)&&(s2 == 2)))
	    {
	      val *= 1.0 - prob_both_on;
	      if(s1 == 2)
		val *= 1.0 - prob_one_off;
	      else
		val *= prob_one_off;
	    }
	  else
	    val = 0.0;

	}


    }
  
  if(val == 0.0)
    {
      cout << "Impossible split probability: " << s<< "->" << s1<<s2 << " p(both on) = "<<prob_both_on << " p_one_off = "<< prob_one_off << endl;
    }
  

  /* else 
    {
      if(s == 1)
	{
	  if((s1 == 1) || (s2 == 1))
	    val = 1.0/5.0;
	}
      else
	{
	  if( ((s1 == 2)&&(s2 != 1)) || ((s2 == 2)&&(s1 != 1)) )
	    val = 1.0/3.0;
	}
	}*/
  return val;
}

void StateMergeSplitUpdater::PrintHyperpars(ostream& str)
{
  short int s, s1, s2;
  str << par_name << " Baseline 2 machines -> 1 machine (0=off,1=spike,2=decay)" << endl;
  str << "  " << setw(6)<< 0 << setw(6) << 1<<setw(6)<< 2  << endl;
  for(s1 = 0; s1 < 3; s1++)
    {
      for(s2 = 0; s2 < 3; s2++)
	{
	  str << s1 << s2 << " ";
	  for(s = 0; s < 3; s++)
	    str << setw(5)<<setprecision(3) << MergeStateProb(s1, s2, s) << " ";
	  str << endl;
	}
    }
}

/*=====================================End State=================================*/



/*=================================Birth/Death/Immune=============================*/


BirthDeathImmuneMergeSplitUpdater::BirthDeathImmuneMergeSplitUpdater(string probfile, vector<double>& hpar1, DistributionList* dists, ParameterUpdater* par):MergeSplitUpdater(hpar1, dists, par)
{
  ifstream infile;
  int i,j;
  double val;
  infile.open(probfile.c_str(), ios::in);
  for(i = 0; i < 36; i++)
    row_sum[i] = 0.0;
  for(i = 0; i < 6; i++)
    col_sum[i] = 0.0;
  for(i = 0; i < 36; i++)
    for(j = 0; j < 6; j++)
	{
	  infile >> val;
	  bdi_prob[i][j] = val;
	  row_sum[i] += val;
	  col_sum[j] += val;
	}
  cur_split_sum = 0;
  cur_merge_sum = 0;
 
  cur_split_prob.clear();
  cur_merge_prob.clear();
}

void BirthDeathImmuneMergeSplitUpdater::PrintHyperpars(ostream& str)
{
  int i;
  short int s1, s2, s;
  short int birth,death,immune;
  int T = 1;
  int last_nonoff=0;
  str << par_name << " Baseline 2 machines -> 1 machine"<<"Format is BDI in (0,1)x(0,1)x(0,1), Birth(0,284), Death(last nonoff,T), Immune(no,yes)"<< endl;
  str << "          000   010   011   100   110   111"<< endl;   
  for(i = 0; i < 36; i++)
    {
      s1=i/6;
      s2=i%6;
      BDIFromIndex(s1,0,1,&birth,&death,&immune);
      birth = (birth > 0);
      str << birth << death<<immune<< ":";
      BDIFromIndex(s2,0,1,&birth,&death,&immune);
      birth = (birth > 0);
      str << birth << death<<immune << " ";
      for(s = 0; s < 6; s++)
	str <<setw(5)<< setprecision(3) << bdi_prob[i][s] << " ";
      str << endl;
    }
}

void BirthDeathImmuneMergeSplitUpdater::MergePropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
{
  //in a merge step, current_machines has m1, m2; and proposed_machines has m
  int bdi1, bdi2, bdi_12;
  int last_non;
  short int newbirth, newdeath, newimmune;
  SetFirstLast(current_machines);
  if(!(CheckFirstLast(proposed_machines)))
    {
      cout << "Error in Merge/Split step:  proposal states do not match current states in first12 = min(first1, first2), last12 = max(last1, last2)" << endl;
    }
  if(last_nonoff[0] > last_nonoff[1])
    last_non = last_nonoff[0];
  else
    last_non = last_nonoff[1];

  bdi1 = current_machines[0]->parameters->BDIindex();
  bdi2 = current_machines[1]->parameters->BDIindex();
  
  SetProbsMerge(bdi1, bdi2);
  bdi_12 = SimpleSample(cur_merge_prob, distributions, cur_merge_sum);
  BDIFromIndex(bdi_12, last_non, proposed_machines[0]->parameters->T, &newbirth, &newdeath, &newimmune);
  //We have sampled a birth/death/immune for the proposed machine (whose states have already been set), now set that value
  proposed_machines[0]->parameters->SetBirthDeathImmune(newbirth, newdeath, newimmune); 
}

double BirthDeathImmuneMergeSplitUpdater::MergeLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
{
  //We need to get the backward compatible probability from the merge step

  double loglike;
  int bdi_12;
  bdi_12 = proposed_machines[0]->parameters->BDIindex();
  SetProbsSplit(bdi_12, first_nonoff[0], first_nonoff[1], 1);
  //cancellation as the p[i][j] index is the same for both row and column
  loglike = log(cur_merge_sum) - log(cur_split_sum);

  return loglike;
}

void BirthDeathImmuneMergeSplitUpdater::SplitPropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
{
  int i;
  int bdi1, bdi2, bdi_12;
  int rowindex;
  short int newbirth1, newdeath1, newimmune1, newbirth2, newdeath2, newimmune2;
  //in a split step, current_machines has m; and proposed_machines has m1, m2
  SetFirstLast(proposed_machines);
  if(!(CheckFirstLast(current_machines)))
    {
      cout << "Error in Merge/Split step:  proposal states do not match current states in first12 = min(first1, first2), last12 = max(last1, last2)" << endl;
    }
  bdi_12 = current_machines[0]->parameters->BDIindex();

  SetProbsSplit(bdi_12, first_nonoff[0], first_nonoff[1]);
  rowindex = SimpleSample(cur_split_prob, distributions, cur_split_sum);
  //cout << "Current probabilities: ";
  //for(i = 0; i < 36; i++)
  //  cout << cur_split_prob[i] << " ";
  //cout << endl << "Current Split sum = " << cur_split_sum << endl;
  //cout << "Sampled row = " << rowindex << endl;
  bdi1 = rowindex/6;
  bdi2 = rowindex%6;
  //cout << "BDI indexes are "<< bdi1 << " " << bdi2 << endl;
  BDIFromIndex(bdi1, last_nonoff[0], proposed_machines[0]->parameters->T, &newbirth1, &newdeath1, &newimmune1);
  BDIFromIndex(bdi2, last_nonoff[1], proposed_machines[1]->parameters->T, &newbirth2, &newdeath2, &newimmune2);
  
  //We have sampled a birth/death/immune for the proposed machines (whose states have already been set), now set those values
  proposed_machines[0]->parameters->SetBirthDeathImmune(newbirth1, newdeath1, newimmune1);
  proposed_machines[1]->parameters->SetBirthDeathImmune(newbirth2, newdeath2, newimmune2);

}

double BirthDeathImmuneMergeSplitUpdater::SplitLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
{
  double loglike;
  int bdi1, bdi2;
  bdi1 = proposed_machines[0]->parameters->BDIindex();
  bdi2 = proposed_machines[1]->parameters->BDIindex();
  SetProbsMerge(bdi1, bdi2, 1);
  //cancellation as the p[i][j] index is the same for both row and column
  loglike = log(cur_split_sum) - log(cur_merge_sum);

  return loglike;
}

void BirthDeathImmuneMergeSplitUpdater::PrintValue(vector<Machine*>& machines, ostream& str)
{
  int i;
  for(i =0; i < machines.size(); i++)
    str << machines[i]->uid << ": ("<<machines[i]->parameters->Birthdate() << " " << machines[i]->parameters->Deathdate() << " " << machines[i]->parameters->Is_immune() << ") ";
  
}


void BirthDeathImmuneMergeSplitUpdater::SetProbsMerge(int bdi1, int bdi2, int sum_only )
{
  int i;
  int rowindex;
  cur_merge_prob.clear();
  if(!sum_only)
    {
      rowindex = bdi1*6 + bdi2;
      for(i = 0; i < 6; i++)
	cur_merge_prob.push_back(bdi_prob[rowindex][i]);
    }
  cur_merge_sum = row_sum[rowindex];
}

void BirthDeathImmuneMergeSplitUpdater::SetProbsSplit(int bdi, int first_nonoff_1, int first_nonoff_2, int sum_only)
{
  int i;
  double cp;
  cur_split_prob.clear();
  //first push back the column needed
  if(!sum_only)
    {
      for(i = 0; i < 36; i++)
	cur_split_prob.push_back(bdi_prob[i][bdi]);
    }
  cur_split_sum = col_sum[bdi];

  //now check the first nonoff state for each of the individuals.  If it does not permit a birth of 284, zero out the prob for the appropriate rows
  //also subtract that probability from the current split sum;
  //Note that it is impossible for bdi > 2 and either first_nonoff_1 < 284 or first_nonoff_2 < 284 
  // because of the way states are merged and split, min(first_1,first_2) >= first_{1+2}
  if(first_nonoff_1 < 284)
    {
      for(i = 18; i < 36; i++)
	{
	  cur_split_sum -= bdi_prob[i][bdi];
	  if(!sum_only)
	    cur_split_prob[i] = 0.0;
	}
    }
  if(first_nonoff_2 < 284)
    {
      {
	for(i = 0; i < 36; i++)
	  {
	    if( (i % 6) > 2)
	      {
		cur_split_sum -= bdi_prob[i][bdi];
		if(!sum_only)
		  cur_split_prob[i] = 0.0;
	      }
	  }
	
      }
    }

  if(cur_split_sum <= 0.0)
    cout << "Error in BirthDeathImmuneMergeSplitUpdater::SetProbsSplit, cur_split_sum = " << cur_split_sum << endl;
}

void BirthDeathImmuneMergeSplitUpdater::SetFirstLast(vector<Machine*>& two_machines)
{
  int i;
  first_nonoff.clear();
  last_nonoff.clear();
  
  for(i = 0; i < 2; i++)
    {
      first_nonoff.push_back(FirstLastRun(two_machines[i]->parameters, 1, 1));
      last_nonoff.push_back(FirstLastRun(two_machines[i]->parameters, -1, 1));
    }
}


int BirthDeathImmuneMergeSplitUpdater::CheckFirstLast(vector<Machine*>& one_machine)
{
  int min_first;
  int max_last;
  int one_first;
  int one_last;
  int valid;
 

  if(first_nonoff[0] < first_nonoff[1])
    min_first = first_nonoff[0];
  else
    min_first = first_nonoff[1];

  if(last_nonoff[0] > last_nonoff[1])
    max_last = last_nonoff[0];
  else
    max_last = last_nonoff[1];
  one_first = FirstLastRun(one_machine[0]->parameters, 1, 1);
  one_last = FirstLastRun(one_machine[0]->parameters, -1, 1);
 
  //cout << "First non-off (1,2,12) = " << setw(4) << first_nonoff[0] << " " << setw(4) <<  first_nonoff[1] <<  " " << setw(4) << one_first << endl; 
  //cout << "Last non-off (1,2,12)  = " << setw(4) << last_nonoff[0] << " " << last_nonoff[1] << " " << setw(4) << one_last << endl; 
  
 
  valid = (one_first == min_first) && (one_last == max_last);
  return(valid);
}

/*=================================End Birth/Death/Immune=============================*/

/*===========================================Q========================================*/

void QBoundedRescaledSampler::gfun1(double q, double &q1, double &q2, vector<double>& proposal_par, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines)
{
  /*
  double A1, A2, A1U2, A1I2;
  if(current_machines.size() == 2)
    ActiveOverlap(current_machines[0]->parameters, current_machines[1]->parameters, A1, A2, A1U2, A1I2);
  else
    ActiveOverlap(proposed_machines[0]->parameters, proposed_machines[1]->parameters, A1, A2, A1U2, A1I2);
  */
  double Ndo, Nod, Ndd;
  double Ydo, Yod, Ydd;
  double pos_d1, pos_d2;
  double pri_mode = 3.0;
  
  Ndo = (double)(y_stats->baselineonlytotal[0]);
  if(Ndo == 0)
    Ndo = 1.0;
  Nod = (double)(y_stats->baselineonlytotal[1]);
  if(Nod == 0)
    Nod = 1.0;
  Ndd = (double)(y_stats->baselineonlytotal[2]);
  if(Ndd == 0)
    Ndd = 1.0;

  Ydo = y_stats->baselineonlymean[0]*Ndo;
  Yod = y_stats->baselineonlymean[1]*Nod;
  Ydd = y_stats->baselineonlymean[2]*Ndd;

  pos_d1 = Ydd/Ndd - Yod/Nod;
  if(pos_d1 < 0)
    pos_d1 = 0;
  pos_d2 = Ydd/Ndd - Ydo/Ndo;
  if(pos_d2 < 0)
    pos_d2 = 0;

  q1 = (1.0/(Ndd + Ndo+ 1.0))*( pri_mode + Ydo + Ndd*pos_d1);
  q2 = (1.0/(Ndd + Nod+ 1.0))*( pri_mode + Yod + Ndd*pos_d1);

  /*
  if(q1 < l)
    {
      q1 = l + 0.001;
    }
  if(q1 > u)
    {
      q1 = u - 0.001;
    }
  if(q2 < l)
    {
      q2 = l + 0.001;
    }
  if(q2 > u)
    {
      q2 = u - 0.001;
    }
  */
}

void QBoundedRescaledSampler::gfun2(double q1, double q2, double& q, vector<double>& proposal_par,vector<Machine*>& current_machines, vector<Machine*>& proposed_machines)
{
  /*double A1, A2, A1U2, A1I2;
  if(current_machines.size() == 2)
    ActiveOverlap(current_machines[0]->parameters, current_machines[1]->parameters, A1, A2, A1U2, A1I2);
  else
    ActiveOverlap(proposed_machines[0]->parameters, proposed_machines[1]->parameters, A1, A2, A1U2, A1I2);
  */
  double Ndecay;
  double Ydecay;
  double pri_mode = 4.4;
  
  Ndecay = (double)(y_stats->baselineonlytotal[3]);
  Ydecay = y_stats->baselineonlymean[3]* Ndecay;
  
  q = 1.0/(Ndecay + 1.0)*(pri_mode + Ydecay);

  /*
  if(q > u)
    {
      q = u - 0.001;
    }
  if(q < l)
    {
      q = l + 0.001;
    }
  */
}

/*
void QMergeSplitUpdater::SetGParameters(vector<Machine*>& current_machines, vector<Machine*>& proposed_machines)
{
  int sz = sampler_pars.size();
  double A1, A2, A1U2, A1I2;
  if(current_machines.size() > 1)
    {
      ActiveOverlap(current_machines[0]->parameters, current_machines[1]->parameters, A1, A2, A1U2, A1I2);
    }
  else
    {
      ActiveOverlap(proposed_machines[0]->parameters, proposed_machines[1]->parameters, A1, A2, A1U2, A1I2);
    }
  y_stats->Initialize(current_machines, proposed_machines);
  sampler_pars[sz-1] = A1;
  sampler_pars[sz-2] = A2;
  sampler_pars[sz-3] = A1U2;
  }*/


QMergeSplitUpdater::QMergeSplitUpdater(vector<double>& hpar1, DistributionList* dists, QUpdater* par):MergeSplitUpdater(hpar1, dists, (ParameterUpdater*)(par))
{
  //l = prior_par[0], u = prior_par[1], prior_a = prior_par[2], prior_b = prior_par[3]
  //prop_alpha = proposal_par[0], samplesize = proposal_par[1],
   double l = prior_hyperpar[0];
   double u = prior_hyperpar[1];
   double samplesize = proposal_par[0];
   int i;
   sampler = new QBoundedRescaledSampler(l, u, samplesize, dists);
 
   /*
  sampler_pars.clear();
   for(i = 1; i < proposal_par.size(); i++)
     sampler_pars.push_back(proposal_par[i]);
  
   sampler_pars.push_back(0.0); //A1
   sampler_pars.push_back(0.0); //A2
   sampler_pars.push_back(0.0); //A1U2
   */
}

void QMergeSplitUpdater::MergePropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
{
 double q1, q2, q12;
  q1 = current_machines[0]->parameters->q;
  q2 = current_machines[1]->parameters->q;
  sampler->y_stats->Initialize(current_machines, proposed_machines);  
  sampler->MergePropose(q1, q2, q12, proposal_par, current_machines, proposed_machines);
  proposed_machines[0]->parameters->q = q12;
}

double QMergeSplitUpdater::MergeLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
{
  double loglike;
  double q1, q2, q12;
  q1 = current_machines[0]->parameters->q;
  q2 = current_machines[1]->parameters->q;
  q12 = proposed_machines[0]->parameters->q;
  loglike = LogPriorRatio(q12, q1, q2, 1) + sampler->MergeLogRatio(q1, q2, q12, proposal_par, current_machines, proposed_machines);
  return loglike;
}

void QMergeSplitUpdater::SplitPropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
{
 double q1, q2, q12;
  q12 = current_machines[0]->parameters->q;
  sampler->y_stats->Initialize(current_machines, proposed_machines);
  sampler->SplitPropose(q12, q1, q2, proposal_par, current_machines, proposed_machines);
  proposed_machines[0]->parameters->q = q1;
  proposed_machines[1]->parameters->q = q2;
}

double QMergeSplitUpdater::SplitLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
{
  double loglike;
  double q1, q2, q12;  
  q1 = proposed_machines[0]->parameters->q;
  q2 = proposed_machines[1]->parameters->q;
  q12 = current_machines[0]->parameters->q;
  loglike = LogPriorRatio(q12, q1, q2, -1) + sampler->SplitLogRatio(q1, q2, q12,proposal_par, current_machines, proposed_machines);
  return loglike;  

}

double QMergeSplitUpdater::LogPriorRatio(double s, double s1, double s2, int mergesplit)
{
  
  double loglike;
  double l = prior_hyperpar[0];
  double u = prior_hyperpar[1];
  double pri_alpha = prior_hyperpar[2];
  double pri_beta =prior_hyperpar[3];
  loglike = LogPriorRatioScaledBeta(s1, s2, s, l, u, pri_alpha, pri_beta);
  return( ((double)(mergesplit))* loglike);
}


void QMergeSplitUpdater::PrintValue(vector<Machine*>& machines, ostream& str)
{
  int i;
  for(i = 0; i < machines.size(); i++)
    {
      str << machines[i]->uid << ": " << machines[i]->parameters->q << " ";
    }
}

/*=======================================End Q========================================*/

/*=========================================Alpha========================================*/
/*
double AlphaBoundedRescaledSampler::gfun1(double q, vector<double>& pars)
{
  return 0.0;
}

double AlphaBoundedRescaledSampler::gfun2(double q1, double q2, vector<double>& pars)
{
  return 0.0;
}

void AlphaMergeSplitUpdater::SetGParameters(vector<Machine*>& current_machines, vector<Machine*>& proposed_machines)
{

}
*/

AlphaMergeSplitUpdater::AlphaMergeSplitUpdater(vector<double>& hpar1, DistributionList* dists, ParameterUpdater* par):MergeSplitUpdater(hpar1, dists, par)
{
  //l = prior_par[0], u = prior_par[1], prior_a = prior_par[2], prior_b = prior_par[3]
   //prop_alpha = proposal_par[0], samplesize = proposal_par[1],
   double l = 0.0;
   double u = 1.0;
   double samplesize = proposal_par[0];
   int i;
   sampler = new AlphaBoundedRescaledSampler(l, u, samplesize, dists);
   /*
   sampler_pars.clear();
   for(i = 1; i < proposal_par.size(); i++)
     sampler_pars.push_back(proposal_par[i]);
   */
   //push back placeholders for any on-the-fly parameter values
}

void AlphaMergeSplitUpdater::MergePropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
{
  double alpha1, alpha2, alpha12;
  alpha1 = current_machines[0]->parameters->alpha;
  alpha2 = current_machines[0]->parameters->alpha;
  sampler->y_stats->Initialize(current_machines, proposed_machines);
  sampler->MergePropose(alpha1, alpha2, alpha12, proposal_par, current_machines, proposed_machines);
  proposed_machines[0]->parameters->alpha = alpha12;
}


double AlphaMergeSplitUpdater::MergeLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
{
  double loglike;
  double alpha1, alpha2, alpha12;
  alpha1 = current_machines[0]->parameters->alpha;
  alpha2 = current_machines[1]->parameters->alpha;
  alpha12 = proposed_machines[0]->parameters->alpha;
  loglike = LogPriorRatio(alpha12, alpha1, alpha2, 1) + sampler->MergeLogRatio(alpha1, alpha2, alpha12, proposal_par, current_machines, proposed_machines);
  return loglike;
}

void AlphaMergeSplitUpdater::SplitPropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
{
  double alpha1, alpha2, alpha12;
  alpha12 = current_machines[0]->parameters->alpha;
  sampler->y_stats->Initialize(current_machines, proposed_machines);
  sampler->SplitPropose(alpha12, alpha1, alpha2, proposal_par, current_machines, proposed_machines);
  proposed_machines[0]->parameters->alpha = alpha1;
  proposed_machines[1]->parameters->alpha = alpha2;
}
  
double AlphaMergeSplitUpdater::SplitLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
{
  double loglike;
  double alpha1, alpha2, alpha12;  
  alpha1 = proposed_machines[0]->parameters->alpha;
  alpha2 = proposed_machines[1]->parameters->alpha;
  alpha12 = current_machines[0]->parameters->alpha;
  loglike = LogPriorRatio(alpha12, alpha1, alpha2, -1) + sampler->SplitLogRatio(alpha1, alpha2, alpha12, proposal_par, current_machines, proposed_machines);
  return loglike;
}
  
double AlphaMergeSplitUpdater::LogPriorRatio(double s, double s1, double s2, int mergesplit)
{
  double loglike;
  double pri_alpha = prior_hyperpar[0];
  double pri_beta =prior_hyperpar[1];
  loglike = LogPriorRatioScaledBeta(s1, s2, s, 0.0, 1.0, pri_alpha, pri_beta);
  return( ((double)(mergesplit))* loglike);
 
}
  
void AlphaMergeSplitUpdater::PrintValue(vector<Machine*>& machines, ostream& str)
{
 int i;
  for(i = 0; i < machines.size(); i++)
    {
      str << machines[i]->uid << ": " << machines[i]->parameters->alpha << " ";
    }
}

/*=======================================End Alpha========================================*/



/*=====================================OmegaScaled========================================*/



void OmegaBoundedRescaledSampler::gfun1(double q, double &q1, double &q2, vector<double>& proposal_par, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines)
{
  double baseline1, baseline2, baseline12;
  double Nso, Nos, Nsd, Nds, Nss;
  double Yso, Yos, Ysd, Yds, Yss;
  double pri_mode = 4.4;
  double all_sum_1;
  double all_sum_2;
  double pos_ss_1, pos_ss_2, pos_sd_1, pos_ds_2;
  Nso = y_stats->spikecounts[0];
  Yso = y_stats->spikesums[0];
  if(Nso == 0)
    Nso = 1.0;

  Nos = y_stats->spikecounts[1];
  Yos = y_stats->spikesums[1];
  if(Nos == 0)
    Nos = 1.0;

  Nsd = y_stats->spikecounts[2];
  Ysd = y_stats->spikesums[2];
  if(Nsd == 0)
    Nsd = 1.0;

  Nds = y_stats->spikecounts[3];
  Yds = y_stats->spikesums[3];
  if(Nds == 0)
    Nds = 1.0;

  Nss = y_stats->spikecounts[4];
  Yss = y_stats->spikesums[4];
  if(Nss == 0)
    Nss = 1.0;

  pos_ss_1 =  (Yss/Nss) - (Yso/Nso);
  if(pos_ss_1 < 0)
    pos_ss_1 = 0;

  pos_ss_2 = (Yss/Nss) - (Yos/Nos);
  if(pos_ss_2 < 0)
    pos_ss_2 = 0;

  pos_sd_1 = (Ysd/Nsd) - baseline2;
  if(pos_sd_1 < 0)
    pos_sd_1 = 0;
  
  pos_ds_2  = (Yds/Nds) - baseline1;
  if(pos_ds_2 < 0)
    pos_ds_2 = 0;

  all_sum_1 = 1.0 + Nso + Nsd + Nss;
  all_sum_2 = 1.0 + Nos + Nds + Nss;

  if(current_machines.size() == 2)
    {
      baseline1 = current_machines[0]->parameters->q;
      baseline2 = current_machines[1]->parameters->q;
      baseline12 = proposed_machines[0]->parameters->q;
    }
  else
    {
      baseline1 = proposed_machines[0]->parameters->q;
      baseline2 = proposed_machines[1]->parameters->q;
      baseline12 = current_machines[0]->parameters->q;
    }

  q1 =  (1.0/all_sum_1)*( pri_mode + Nss*( pos_ss_1/baseline1)  + Nsd*( pos_sd_1/baseline1) + Nso/baseline1 );
  q2 =  (1.0/all_sum_2)*( pri_mode + Nss*( pos_ss_2/baseline2)  + Nds*( pos_ds_2/baseline2) + Nos/baseline2 ); 
  
}

void OmegaBoundedRescaledSampler::gfun2(double q1, double q2, double& q, vector<double>& proposal_par,vector<Machine*>& current_machines, vector<Machine*>& proposed_machines)
{
  double baseline12;
  double pri_mode = 4.4;
  if(current_machines.size() == 2)
    {
      baseline12 = proposed_machines[0]->parameters->q;
    }
  else
    {
      baseline12 = current_machines[0]->parameters->q;
    }

  q = (1.0/((double)(y_stats->spikecounttotal[3])+1.0)) *(pri_mode +  (double)(y_stats->spikecounttotal[3])* y_stats->spikecountmean[3]/(baseline12));
 
}




OmegaScaledMergeSplitUpdater::OmegaScaledMergeSplitUpdater(vector<double>& hpar1, DistributionList* dists, ParameterUpdater* par):MergeSplitUpdater(hpar1, dists, par)
 {
   //l = prior_par[0], u = prior_par[1], prior_a = prior_par[2], prior_b = prior_par[3]
   //prop_alpha = proposal_par[0], samplesize = proposal_par[1],
    l = prior_hyperpar[0];
    u = prior_hyperpar[1];
    samplesize = proposal_par[0];
   int i;

   sampler = new OmegaBoundedRescaledSampler(l, u, samplesize, dists);
   
   //sampler_pars.clear();
   //for(i = 1; i < proposal_par.size(); i++)
   //  sampler_pars.push_back(proposal_par[i]);
   
   //push back placeholders for any on-the-fly parameter values
 }

void OmegaScaledMergeSplitUpdater::MergePropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
{
  double omega1, omega2, omega12;
  double min_u, max_u;
  double baseline1, baseline2, baseline12;
  double baseline_ssd, baseline_remain;

  omega1 = current_machines[0]->parameters->omega;
  omega2 = current_machines[1]->parameters->omega;
  sampler->y_stats->Initialize(current_machines, proposed_machines);

  

  sampler->MergePropose(omega1, omega2, omega12, proposal_par, current_machines, proposed_machines);



  proposed_machines[0]->parameters->omega = omega12;
}
 
double OmegaScaledMergeSplitUpdater::MergeLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
{
  double loglike;
  double omega1, omega2, omega12;
  omega1 = current_machines[0]->parameters->omega;
  omega2 = current_machines[1]->parameters->omega;
  omega12 = proposed_machines[0]->parameters->omega;
  loglike = LogPriorRatio(omega12, omega1, omega2, 1) + sampler->MergeLogRatio(omega1, omega2, omega12,proposal_par, current_machines, proposed_machines);
  return loglike;
}
 

void OmegaScaledMergeSplitUpdater::SplitPropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
 {
   double omega1, omega2, omega12;
   double min_u, max_u;
   double baseline1, baseline2, baseline12;
   double baseline_ssd, baseline_remain;
   
   omega12 = current_machines[0]->parameters->omega;
 
   sampler->y_stats->Initialize(current_machines, proposed_machines);
   
   sampler->SplitPropose(omega12, omega1, omega2, proposal_par, current_machines, proposed_machines);
   proposed_machines[0]->parameters->omega = omega1;
   proposed_machines[1]->parameters->omega = omega2;
 }

double OmegaScaledMergeSplitUpdater::SplitLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
{
  double loglike;
  double omega1, omega2, omega12;  
  omega1 = proposed_machines[0]->parameters->omega;
  omega2 = proposed_machines[1]->parameters->omega;
  omega12 = current_machines[0]->parameters->omega;
  loglike = LogPriorRatio(omega12, omega1, omega2, -1) + sampler->SplitLogRatio(omega1, omega2, omega12, proposal_par, current_machines, proposed_machines);
  return loglike;
}

double OmegaScaledMergeSplitUpdater::LogPriorRatio(double s, double s1, double s2, int mergesplit)
{
  double loglike;
  double l = prior_hyperpar[0];
  double u = prior_hyperpar[1];
  double pri_alpha = prior_hyperpar[2];
  double pri_beta =prior_hyperpar[3];
  loglike = LogPriorRatioScaledBeta(s1, s2, s, l, u, pri_alpha, pri_beta);
  return( ((double)(mergesplit))* loglike);
}

void OmegaScaledMergeSplitUpdater::PrintValue(vector<Machine*>& machines, ostream& str)
{
  int i;
  for(i = 0; i < machines.size(); i++)
    str << machines[i]->uid << ": " << machines[i]->parameters->omega << " ";
}

OmegaScaledMergeSplitUpdater::~OmegaScaledMergeSplitUpdater()
{
  delete sampler;
}


/*==================================End OmegaScaled========================================*/


/*========================================Offlambda========================================*/

void OffLambdaMergeSplitUpdater::MergePropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
{
  double merge_k_val;
  double merge_theta_val;
  double prop_offlambda;
  double pri_k = prior_hyperpar[0];
  double pri_theta = prior_hyperpar[1];
  merge_k_val = (double)(proposed_machines[0]->parameters->NumOffOffTrans()) + pri_k;
  merge_theta_val = 1.0/((double)(proposed_machines[0]->parameters->NumOffSeq()) +  (1.0/pri_theta));
  distributions->Simulate(4,1, &prop_offlambda, merge_k_val, merge_theta_val);
  proposed_machines[0]->parameters->off_lambda = prop_offlambda;

}

double OffLambdaMergeSplitUpdater::MergeLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
{
  double merge_k_val;
  double merge_theta_val;
  double prop_offlambda;
  double s1_k_val, s2_k_val;
  double s1_theta_val, s2_theta_val;
  double s1_offlambda, s2_offlambda;
  double logratio;
  double pri_k = prior_hyperpar[0];
  double pri_theta = prior_hyperpar[1];
  
  s1_offlambda = current_machines[0]->parameters->off_lambda;
  s2_offlambda = current_machines[1]->parameters->off_lambda;
  prop_offlambda = proposed_machines[0]->parameters->off_lambda;
  logratio = LogPriorRatio(s1_offlambda, s2_offlambda, prop_offlambda, 1);

  merge_k_val = (double)(proposed_machines[0]->parameters->NumOffOffTrans()) + pri_k;
  merge_theta_val = 1.0/((double)(proposed_machines[0]->parameters->NumOffSeq()) +  (1.0/pri_theta));
  s1_k_val = (double)(current_machines[0]->parameters->NumOffOffTrans()) + pri_k;
  s1_theta_val = 1.0/((double)(current_machines[0]->parameters->NumOffSeq()) +  (1.0/pri_theta));
  s2_k_val = (double)(current_machines[1]->parameters->NumOffOffTrans()) + pri_k;
  s2_theta_val = 1.0/((double)(current_machines[1]->parameters->NumOffSeq()) +  (1.0/pri_theta));

  logratio += distributions->LogDensity(4, s1_offlambda, s1_k_val, s1_theta_val) + distributions->LogDensity(4, s2_offlambda, s2_k_val, s2_theta_val) - distributions->LogDensity(4, prop_offlambda, merge_k_val, merge_theta_val);

  return logratio;
}

void OffLambdaMergeSplitUpdater::SplitPropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
{
  double s1_k_val, s2_k_val;
  double s1_theta_val, s2_theta_val;
  double s1_offlambda, s2_offlambda;
  double pri_k = prior_hyperpar[0];
  double pri_theta = prior_hyperpar[1];
  // cout << "Prior K = " << pri_k << "Prior Theta = " << pri_theta << endl;
  s1_k_val = (double)(proposed_machines[0]->parameters->NumOffOffTrans()) + pri_k;
  s1_theta_val = 1.0/((double)(proposed_machines[0]->parameters->NumOffSeq()) +  (1.0/pri_theta));
  distributions->Simulate(4,1, &s1_offlambda, s1_k_val, s1_theta_val);
  proposed_machines[0]->parameters->off_lambda = s1_offlambda;

  s2_k_val = (double)(proposed_machines[1]->parameters->NumOffOffTrans()) + pri_k;
  s2_theta_val = 1.0/((double)(proposed_machines[1]->parameters->NumOffSeq()) +  (1.0/pri_theta));
  distributions->Simulate(4,1, &s2_offlambda, s2_k_val, s2_theta_val);
  proposed_machines[1]->parameters->off_lambda = s2_offlambda;

}
 
double OffLambdaMergeSplitUpdater::SplitLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
{
  double merge_k_val;
  double merge_theta_val;
  double prop_offlambda;
  double s1_k_val, s2_k_val;
  double s1_theta_val, s2_theta_val;
  double s1_offlambda, s2_offlambda;
  double logratio;
  double pri_k = prior_hyperpar[0];
  double pri_theta = prior_hyperpar[1];

  s1_offlambda = proposed_machines[0]->parameters->off_lambda;
  s2_offlambda = proposed_machines[1]->parameters->off_lambda;
  prop_offlambda = current_machines[0]->parameters->off_lambda;
  logratio = LogPriorRatio(s1_offlambda, s2_offlambda, prop_offlambda, -1);
  //cout << "LogPriorRatio for offlambda: " << logratio<<endl;
  merge_k_val = (double)(current_machines[0]->parameters->NumOffOffTrans()) + pri_k;
  merge_theta_val = 1.0/((double)(current_machines[0]->parameters->NumOffSeq()) +  (1.0/pri_theta));
  s1_k_val = (double)(proposed_machines[0]->parameters->NumOffOffTrans()) + pri_k;
  s1_theta_val = 1.0/((double)(proposed_machines[0]->parameters->NumOffSeq()) +  (1.0/pri_theta));
  s2_k_val = (double)(proposed_machines[1]->parameters->NumOffOffTrans()) + pri_k;
  s2_theta_val = 1.0/((double)(proposed_machines[1]->parameters->NumOffSeq()) +  (1.0/pri_theta));

  logratio += distributions->LogDensity(4, prop_offlambda, merge_k_val, merge_theta_val) - distributions->LogDensity(4, s1_offlambda, s1_k_val, s1_theta_val) - distributions->LogDensity(4, s2_offlambda, s2_k_val, s2_theta_val) ;
  //cout << "posterior log ratio : " << logratio << endl;

  return logratio;

}

double OffLambdaMergeSplitUpdater::LogPriorRatio(double s1, double s2, double s, int mergesplit)
{
  //default is log(12) - log(1,2) eg a mergestep 
  //mergesplit = 1: merge, mergesplit = -1, split
  //a = v1, b = v2, c = v
  double loglike;
  loglike = LogPriorRatioScaledGamma(s1, s2, s, 0, -1, prior_hyperpar[0], prior_hyperpar[1]);
  return((double)(mergesplit)*loglike);
}
  
void OffLambdaMergeSplitUpdater::PrintValue(vector<Machine*>& machines, ostream& str)
{
  int i;
  for(i = 0; i < machines.size(); i++)
    {
      str << machines[i]->uid << ": " << machines[i]->parameters->off_lambda << " ";
    }
}
  


/*=======================================End Offlambda==================================*/



/*===========================================Rho========================================*/

RhoMergeSplitUpdater::RhoMergeSplitUpdater( vector<double>& hpar, DistributionList* dists, RhoUpdater* par):MergeSplitUpdater(hpar, dists, (ParameterUpdater*)(par))
{
  state = par->state;
}

void RhoMergeSplitUpdater::MergePropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
{
  //Assumes nu has already been set
  double rho12;
  rho12_mean = LeastSquares(proposed_machines[0], rho12_ss);
  distributions->Simulate(2, 1, &rho12, rho12_mean*rho12_ss, (1.0-rho12_mean)*rho12_ss);
  proposed_machines[0]->parameters->rho[state] = rho12;
  
}

double RhoMergeSplitUpdater::MergeLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
{
  double loglike;
  double rho1, rho2, rho12;
  rho1_mean = LeastSquares(current_machines[0], rho1_ss);
  rho2_mean = LeastSquares(current_machines[1], rho2_ss);
  rho1 = current_machines[0]->parameters->rho[state];
  rho2 = current_machines[1]->parameters->rho[state];
  rho12= proposed_machines[0]->parameters->rho[state];

  //12 - 1 - 2
  //cout << "Rho value(1,2,12) = " << rho1 << " " << rho2 << " " << rho12 << endl;
  //cout << "Rho_mean (1,2,12) = " << rho1_mean << " " << rho2_mean << " " << rho12_mean << endl;
  //cout << "Rho_ss   (1,2,12) = " << rho1_ss << " "<< rho2_ss << " " << rho12_ss << endl;
  //cout << "Rho_alpha(1,2,12) = " << rho1_mean*rho1_ss << " " << rho2_mean * rho2_ss << " " << rho12_mean*rho12_ss << endl;
  //cout << "Rho_beta (1,2,12) = " << (1-rho1_mean)*rho1_ss << " "<< (1-rho2_mean)* rho2_ss << " " << (1-rho12_mean)*rho12_ss << endl;

  loglike = distributions->LogDensity(2, rho1, rho1_mean*rho1_ss, (1.0-rho1_mean)*rho1_ss) + distributions->LogDensity(2, rho2, rho2_mean*rho2_ss, (1.0-rho2_mean)*rho2_ss) - distributions->LogDensity(2, rho12, rho12_mean*rho12_ss, (1.0-rho12_mean)*rho12_ss);


  return loglike;
}
 
void RhoMergeSplitUpdater::SplitPropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
{
  double rho1, rho2;
  rho1_mean = LeastSquares(proposed_machines[0], rho1_ss);
  distributions->Simulate(2, 1, &rho1, rho1_mean*rho1_ss, (1.0-rho1_mean)*rho1_ss);
  proposed_machines[0]->parameters->rho[state] = rho1;

  rho2_mean = LeastSquares(proposed_machines[1], rho2_ss);
  distributions->Simulate(2, 1, &rho2, rho2_mean*rho2_ss, (1.0-rho2_mean)*rho2_ss);
  proposed_machines[1]->parameters->rho[state] = rho2;  
}

double RhoMergeSplitUpdater::SplitLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
{
  double loglike;
  double rho1, rho2, rho12;
  rho12_mean = LeastSquares(current_machines[0], rho12_ss);
  rho12 = current_machines[0]->parameters->rho[state];
  rho1 = proposed_machines[0]->parameters->rho[state];
  rho2 = proposed_machines[1]->parameters->rho[state];
  //12 - 1 - 2
  //cout << "Rho value(1,2,12) = " << rho1 << " " << rho2 << " " << rho12 << endl;
  //cout << "Rho_mean (1,2,12) = " << rho1_mean << " " << rho2_mean << " " << rho12_mean << endl;
  //cout << "Rho_ss   (1,2,12) = " << rho1_ss << " "<< rho2_ss << " " << rho12_ss << endl;
  //cout << "Rho_alpha(1,2,12) = " << rho1_mean*rho1_ss << " " << rho2_mean * rho2_ss << " " << rho12_mean*rho12_ss << endl;
  //cout << "Rho_beta (1,2,12) = " << (1-rho1_mean)*rho1_ss << " "<< (1-rho2_mean)* rho2_ss << " " << (1-rho12_mean)*rho12_ss << endl;
  loglike = distributions->LogDensity(2, rho12, rho12_mean*rho12_ss, (1.0-rho12_mean)*rho12_ss) - distributions->LogDensity(2, rho1, rho1_mean*rho1_ss, (1.0-rho1_mean)*rho1_ss) - distributions->LogDensity(2, rho2, rho2_mean*rho2_ss, (1.0-rho2_mean)*rho2_ss);

  return loglike;
}

double RhoMergeSplitUpdater::LogPriorRatio(double s1, double s2, double s, int mergesplit)
{
   //default is log(12) - log(1,2) eg a mergestep 
  //mergesplit = 1: merge, mergesplit = -1, split
  //a = v1, b = v2, c = v
  double loglike;
  double pri_alpha = prior_hyperpar[0];
  double pri_beta = prior_hyperpar[1];
  loglike = LogPriorRatioScaledBeta(s1, s2, s, 0.0, 1.0, pri_alpha, pri_beta);
  return((double)(mergesplit)*loglike);

}
 
void RhoMergeSplitUpdater::PrintValue(vector<Machine*>& machines, ostream& str)
{
 int i;
  for(i = 0; i < machines.size(); i++)
    {
      str << machines[i]->uid << ": " << machines[i]->parameters->rho[state] << " ";
    }
}

double RhoMergeSplitUpdater::LeastSquares(Machine* machine, double& samplesize)
{
  int successes, totals;
  double rho_hat;
  double a_i;
  double y_i;
  double ssq_i;
  double nu;
  double tstar;
  double num, denom;
  int i;
  double pri_mean;
  double pri_alpha = prior_hyperpar[0];
  double pri_beta = prior_hyperpar[1];
  nu = machine->parameters->nu[state];
  tstar = (double)(machine->parameters->tstar);
  pri_mean = pri_alpha/(pri_alpha + pri_beta);
  samplesize = (pri_alpha + pri_beta);
  RhoNuSuccessesAndTotalCounts(state, machine->parameters, 5, successes, totals);
  rho_hat = (1.0/( (double)(totals) + 1.0)) *( pri_mean + (double)(successes));
  /*
  num = 0.0;
  denom = 0.0;
  for(i = 0; i < 24; i++)
    {
      RhoNuSuccessesAndTotalCounts(state, machine->parameters, i, successes, totals);
      a_i = (sin(2.0*M_PI*(i - tstar)/24) + nu + 1.0)/(nu + 2.0);
      y_i = ((1.0)/((double)(totals)+((pri_alpha+pri_beta)))) * ( pri_mean*a_i  +  ((double)(successes)));
      if(y_i < MIN_LOG_P)
	{
	  y_i = MIN_LOG_P;
	}
      ssq_i = (y_i)*(1.0-y_i)/( (double)(totals)+(pri_alpha+pri_beta) );
      num += (a_i*y_i)/ssq_i;
      denom += (a_i/ssq_i);
      samplesize += (double)(totals);
    }
  if(denom < MIN_LOG_P)
    {
      denom = MIN_LOG_P;
    }
  if(num < MIN_LOG_P)
    {
      num= MIN_LOG_P;
    }
  rho_hat = exp(log(num) - log(denom));
  */
  if(rho_hat < MIN_LOG_P)
    rho_hat = MIN_LOG_P;
  if(rho_hat > 1.0 - MIN_LOG_P)
    rho_hat = 1.0 - MIN_LOG_P;
  return(rho_hat);
}


/*=======================================End Rho========================================*/



/*========================================Gamma========================================*/

 
void GammaMergeSplitUpdater::MergePropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
{
  double prop_gamma;
  double pri_alpha = prior_hyperpar[0];
  double pri_beta = prior_hyperpar[1];
  GetPlusMinusState(state, proposed_machines[0]->parameters, plusstatecount_12, minusstatecount_12);
  distributions->Simulate(2, 1, &prop_gamma, pri_alpha + plusstatecount_12, pri_beta + minusstatecount_12);
  prop_gamma = BetaConstrain(prop_gamma);
  proposed_machines[0]->parameters->gamma[state] = prop_gamma;

}

double GammaMergeSplitUpdater::MergeLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
{
  double gamma_s1, gamma_s2, gamma_s12;
  double loglike;
  double pri_alpha = prior_hyperpar[0];
  double pri_beta = prior_hyperpar[1];

  GetPlusMinusState(state, current_machines[0]->parameters, plusstatecount_1, minusstatecount_1);
  GetPlusMinusState(state, current_machines[1]->parameters, plusstatecount_2, minusstatecount_2);
  gamma_s1 = current_machines[0]->parameters->gamma[state];
  gamma_s2 = current_machines[1]->parameters->gamma[state];
  gamma_s12 = proposed_machines[0]->parameters->gamma[state];
  
  loglike = LogPriorRatio(gamma_s1, gamma_s2, gamma_s12, 1);
  loglike += distributions->LogDensity(2, gamma_s1,pri_alpha+plusstatecount_1, pri_beta+minusstatecount_1)
    +distributions->LogDensity(2, gamma_s2, pri_alpha+plusstatecount_2, pri_beta + minusstatecount_1)
    -distributions->LogDensity(2, gamma_s12, pri_alpha+plusstatecount_12, pri_beta + minusstatecount_12);

  return loglike;
}

void GammaMergeSplitUpdater::SplitPropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
{
  double s1_gamma, s2_gamma;
  double pri_alpha = prior_hyperpar[0];
  double pri_beta = prior_hyperpar[1];

  GetPlusMinusState(state, proposed_machines[0]->parameters, plusstatecount_1, minusstatecount_1);
  GetPlusMinusState(state, proposed_machines[1]->parameters, plusstatecount_2, minusstatecount_2);
  distributions->Simulate(2,1,&s1_gamma, pri_alpha + plusstatecount_1, pri_beta + minusstatecount_1);
  distributions->Simulate(2,1,&s2_gamma, pri_alpha + plusstatecount_2, pri_beta + minusstatecount_2);
  s1_gamma = BetaConstrain(s1_gamma);
  s2_gamma = BetaConstrain(s2_gamma);
  proposed_machines[0]->parameters->gamma[state] = s1_gamma;
  proposed_machines[1]->parameters->gamma[state] = s2_gamma;

}
  
double GammaMergeSplitUpdater::SplitLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
{
  double gamma_s1, gamma_s2, gamma_s12;
  double loglike;
  double pri_alpha = prior_hyperpar[0];
  double pri_beta = prior_hyperpar[1];
  
  GetPlusMinusState(state, current_machines[0]->parameters, plusstatecount_12, minusstatecount_12);
  gamma_s1 = proposed_machines[0]->parameters->gamma[state];
  gamma_s2 = proposed_machines[1]->parameters->gamma[state];
  gamma_s12 = current_machines[0]->parameters->gamma[state];
  
  loglike = LogPriorRatio(gamma_s1, gamma_s2, gamma_s12, -1);
  loglike += distributions->LogDensity(2, gamma_s12, pri_alpha+plusstatecount_12, pri_beta + minusstatecount_12)  - distributions->LogDensity(2, gamma_s1,pri_alpha +plusstatecount_1, pri_beta +minusstatecount_1) - distributions->LogDensity(2, gamma_s2, pri_alpha +plusstatecount_2, pri_beta + minusstatecount_1);

  return loglike;
}
  
double GammaMergeSplitUpdater::LogPriorRatio(double s1, double s2, double s, int mergesplit)
{
  //default is log(12) - log(1,2) eg a mergestep 
  //mergesplit = 1: merge, mergesplit = -1, split
  //a = v1, b = v2, c = v
  double loglike;
  double pri_alpha = prior_hyperpar[0];
  double pri_beta = prior_hyperpar[1];
  
  loglike = LogPriorRatioScaledBeta(s1, s2, s, 0, 1, pri_alpha, pri_beta);
  return((double)(mergesplit)*loglike);
}

void GammaMergeSplitUpdater::PrintValue(vector<Machine*>& machines, ostream& str)
{
  int i;
  for(i = 0; i < machines.size(); i++)
    {
      str << machines[i]->uid << ": " << machines[i]->parameters->gamma[state] << " "; 
    }
}


/*=======================================End Gamma========================================*/



/*===========================================Nu========================================*/

NuMergeSplitUpdater::NuMergeSplitUpdater( vector<double>& hpar, DistributionList* dists, NuUpdater* par):MergeSplitUpdater(hpar, dists, (ParameterUpdater*)(par))
{
  state = par->state;
  double l, u, alpha, beta;
  double sampsize;
  l = prior_hyperpar[0];
  u = prior_hyperpar[1];
  alpha = prior_hyperpar[2];
  beta = prior_hyperpar[3];
  sampsize = proposal_par[0];
  sampler_pars.clear();
  sampler= new BoundedRescaledSampler(l, u, sampsize, dists );
}
 

void NuMergeSplitUpdater::MergePropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
  {
    double alpha, beta;
    double nu1, nu2, nu12;
    nu1 = current_machines[0]->parameters->nu[state];
    nu2 = current_machines[1]->parameters->nu[state];    
    sampler->MergePropose(nu1, nu2, nu12, proposal_par, current_machines, proposed_machines);
    proposed_machines[0]->parameters->nu[state] = nu12;

  }

double NuMergeSplitUpdater::MergeLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
 {
   double loglike;
  
   double nu1, nu2, nu12;
   nu1 = current_machines[0]->parameters->nu[state];
   nu2 = current_machines[1]->parameters->nu[state];
   nu12 = proposed_machines[0]->parameters->nu[state];
  
   loglike = LogPriorRatio(nu1, nu2, nu12, 1) + sampler->MergeLogRatio(nu1, nu2, nu12, proposal_par, current_machines, proposed_machines);
   
   return loglike;
 }

void NuMergeSplitUpdater::SplitPropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
{
  double alpha, beta;
  double nu1, nu2, nu12;
  nu12 = current_machines[0]->parameters->nu[state];
  sampler->SplitPropose(nu12, nu1, nu2, proposal_par, current_machines, proposed_machines);
  proposed_machines[0]->parameters->nu[state] = nu1;
  proposed_machines[1]->parameters->nu[state] = nu2;
}

double NuMergeSplitUpdater::SplitLogRatio(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
{

  double loglike;
  
   double nu1, nu2, nu12;
   nu1 = proposed_machines[0]->parameters->nu[state];
   nu2 = proposed_machines[1]->parameters->nu[state];
   nu12 = current_machines[0]->parameters->nu[state];

   loglike = LogPriorRatio(nu1, nu2, nu12, -1) + sampler->SplitLogRatio(nu1, nu2, nu12, proposal_par, current_machines, proposed_machines);
   
   return loglike;
}

double NuMergeSplitUpdater::LogPriorRatio(double s1, double s2, double s, int mergesplit)
{
  double loglike;
  double l = prior_hyperpar[0];
  double u = prior_hyperpar[1];
  double pri_alpha = prior_hyperpar[2];
  double pri_beta = prior_hyperpar[3];
  loglike = LogPriorRatioScaledBeta(s1, s2, s,l, u, pri_alpha, pri_beta );
  return((double)(mergesplit)*loglike);
}

void NuMergeSplitUpdater::PrintValue(vector<Machine*>& machines, ostream& str)
 {
  int i;
  for(i = 0; i < machines.size(); i++)
    {
      str << machines[i]->uid << ": " << machines[i]->parameters->nu[state] << " "; 
    }
 }

/*=======================================End Nu========================================*/


/*========================================Info========================================*/


void InfoMergeSplitUpdater::MergePropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
{
  int tstar1, tstar2, tstar12;
  tstar1 = current_machines[0]->data->tstar;
  tstar2 = current_machines[1]->data->tstar;
  tstar12 = proposed_machines[0]->data->tstar;
  
  if( tstar1 != tstar2)
    cout << "Unmatching T-stars in merge components (shouldn't happen yet)" << endl;
  if( (tstar12 != tstar1) || (tstar12 != tstar2))
    cout << "Merged T-star not equal to current t-stars: " << tstar12 << " " << tstar1 << " " << tstar2;

}

void InfoMergeSplitUpdater::SplitPropose(vector<Network*>& networks, vector<Machine*>& current_machines, vector<Machine*>& proposed_machines, global_par_obj* global_pars)
{
  int tstar1, tstar2, tstar12;
  tstar1 = proposed_machines[0]->data->tstar;
  tstar2 = proposed_machines[1]->data->tstar;
  tstar12 = current_machines[0]->data->tstar;
  
  if( tstar1 != tstar2)
    cout << "Unmatching T-stars in splitted components (shouldn't happen yet)" << endl;
  if( (tstar12 != tstar1) || (tstar12 != tstar2))
    cout << "current T-star not equal to split t-stars: " << tstar12 << " " << tstar1 << " " << tstar2;

}

void InfoMergeSplitUpdater::PrintValue(vector<Machine*>& machines, ostream& str)
{
  int i;
  for(i = 0; i < machines.size(); i++)
    str << machines[i]->uid << ": (" << machines[i]->data->T << " " << machines[i]->data->timezone << " " << machines[i]->data->tstar << " " << machines[i]->data->iteration_born<< ") ";

}
  


/*=======================================End Info========================================*/


/* Counts are in countupdaters.h */

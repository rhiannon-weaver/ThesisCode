#include "../networkobjects.h"

int ip_dist(unsigned int a, unsigned int b)
{
  int d = 0;
  while(a != b)
    {
      d++;
      a >> 1;
      b >> 1;
    }
  return(d);
}

bool Network::operator<(Network rhs)
{
  return(data->block_id < (rhs.data)->block_id);
}

Network::Network(int N, int rid, netw_par_obj* netw_pars)
{
  data = new netw_dat_obj(N, rid);
  parameters = new netw_par_obj(netw_pars);
}

Network::Network(int N, int rid, ifstream &spikeinfo , ifstream &flinfo, ifstream &ipinfo, netw_par_obj* netw_pars)
{
  data = new netw_dat_obj(N, rid, spikeinfo, flinfo, ipinfo);
  parameters = new netw_par_obj(netw_pars);
}

Network::Network(int N, int rid, ifstream &spikeinfo, ifstream &flinfo, ifstream &ipinfo)
{
  data = new netw_dat_obj(N, rid, spikeinfo, flinfo, ipinfo);
  parameters = new netw_par_obj();
}

Network::Network(Network* cpy)
{
  data = new netw_dat_obj(cpy->data);
  parameters = new netw_par_obj(cpy->parameters);
}

void Network::Print(ostream& str, int N)
{
  str << "%%%%%%%%%%%%%%%  Network " << data->block_id << " %%%%%%%%%%%%%%%%" << endl;
  parameters->Print(str);
  data->Print(str,N);
  cout << "%%%%%%%%%%%%%%%  End Network " << data->block_id << " %%%%%%%%%%%%%%%%"<< endl;
}

Network::~Network()
{
  delete data;
  delete parameters;
}

machine_id_t Machine::id = 0;

machine_id_t Machine::getID()
{
  return id++;
}

void Machine::ResetIDs()
{
  id = 0;
}

Machine::Machine(Network* initial_network,  double* psi_params, double q_start, double alpha_start, double omega_start, double* rho_start, double* gamma_start, double * nu_start, double off_min, double off_max, int immune)
{
  int converge;
  uid = getID();
  data = new mach_dat_obj(initial_network->data, uid, 0);
  parameters = new curr_par_obj(data, psi_params, q_start, alpha_start, omega_start, rho_start, gamma_start, nu_start, off_min, off_max, immune);
  par_acc_allocated = 0;
 
}

Machine::Machine(Network* initial_network, curr_par_obj* initial_pars, hype_par_obj* hyperpars)
{
  int converge;
 uid = getID();
 data = new mach_dat_obj(initial_network->data, uid, 0);
 parameters = new curr_par_obj(data, initial_pars, hyperpars);
 par_acc_allocated = 0;
}

int Machine::RefineStartValuesUsingEM(ostream& printstream)
{
  int converge;
  vector<int> fl;
  vector<short int> st;
  int i;
  for(i = 0; i < data->T; i++)
    {
      fl.push_back(data->GetTotalCount(i));
      st.push_back(parameters->State(i));
    }
  converge = DataDrivenStartValRefine(&(fl), &(st), &(parameters->q), &(parameters->alpha), &(parameters->omega), parameters->rho, parameters->gamma, parameters->nu, printstream);
  for(i = 0; i < parameters->T; i++)
    {
      parameters->SetState(i, st[i]);
    }
}


Machine::Machine(Machine* cpy, int iteration, int get_new_id)
{
  int i;
  data= new mach_dat_obj(cpy->data);
  parameters = new curr_par_obj(cpy->parameters);
  if(get_new_id)
    {
      uid = getID();
      data->uid = uid;
      data->iteration_born = iteration;
      par_acc_allocated = 0;
      num_par = 0;
    }
  else
    {
      uid = cpy->uid;
      par_acc_allocated = cpy->par_acc_allocated;
      num_par = cpy->num_par;
      if(par_acc_allocated)
	{
	  par_acc = new int[num_par];
	  par_total = new int[num_par];
	  for(i = 0; i < num_par; i++)
	    {
	      par_acc[i] = cpy->par_acc[i];
	      par_total[i] = cpy->par_total[i];
	    }
	}
    }
}

void Machine::Print(ostream& str, int N)
{
  str << "%%%%%%%%%%%%%%%  Machine " << uid << " %%%%%%%%%%%%%%%%" << endl;
  data->Print(str, N);
  parameters->Print(str, N);
  cout << "%%%%%%%%%%%%%%%  End Machine " << uid << " %%%%%%%%%%%%%%%%" << endl;
}

double Machine::GetRate(int t)
{
  return(parameters->GetRate(t));
}

Machine::~Machine()
{
  delete data;
  delete parameters;				
  if(par_acc_allocated)
    {
      delete par_acc;
      delete par_total;
    }
		      
}


YTotalActiveStats::YTotalActiveStats(int bl )
{
  y_min = -1;
  y_max = -1;
  y_num_active = 0;
  baseline_lag = bl;
  density.clear();
  cum_sum.clear();
}

void YTotalActiveStats::Initialize( vector<Machine*>& machines)
{
  int i,j;
  int cur_y;
  y_min = -1;
  y_max = -1;
  y_num_active = 0;
  vector<int> counts;
  density.clear();
  cum_sum.clear();
  int is_active;
  
  for(i = 0; i < machines[0]->parameters->T; i++)
    {
      cur_y = machines[0]->data->GetTotalCount(i);
      is_active = machines[0]->parameters->State(i) != 0;
      for(j = 1; j < machines.size(); j++)
	{
	  cur_y += machines[j]->data->GetTotalCount(i);
	  is_active = machines[j]->parameters->State(i) != 0;
	}

      if(is_active)
	{
	  y_num_active++;
	  if(y_min == -1)
	    {
	      y_min = cur_y;
	    }
	  else
	    {
	      if(cur_y < y_min)
		y_min = cur_y;
	    }
	  if(y_max == -1)
	    {
	      y_max = cur_y;
	    }
	  else
	    {
	      if(cur_y > y_max)
		y_max = cur_y;
	    }
	  counts.push_back(cur_y);
	}
    }
  for(i = y_min; i <= y_max; i++)
    {
      density.push_back(0);
      cum_sum.push_back(0);
    }
  for(j = 0; j < counts.size(); j++)
    {
      density[counts[j] - y_min]++;
    }
  cum_sum[0] = density[0];
  for(j = 1; j < cum_sum.size(); j++)
    {
      cum_sum[j] = cum_sum[j-1] + density[j];
    }
  counts.clear();
}

double YTotalActiveStats::YCDF(int y)
{
  double val;
  if( y < y_min)
    {
      val = 0.0;
    }
  else
    {
      if(y > y_max)
	{
	  val = 1.0;
	}
      else
	{
	  val = (double)(cum_sum[y - y_min])/(1.01*(double)(y_num_active));
	}
    }
  return val;
}
 
double YTotalActiveStats::YDensity(int y)
{
  double val;
  if( y < y_min)
    {
      val = 0.0;
    }
  else
    {
      if(y > y_max)
	{
	  val = 0.0;
	}
      else
	{
	  val = (double)(density[y - y_min])/(double)(y_num_active);
	}
    }
  return val;
}
  
/*
void YStateAndCountStats::SetPrior(vector<Machine*>& current_machines)
{
  double pri_q, pri_omega;

  if(current_machines.size() == 1)
    {
      pri_q = current_machines[0]->parameters->q;
      pri_omega = current_machines[0]->parameters->omega;
    }
  else
    {
      pri_q = current_machines[0]->parameters->q + current_machines[1]->parameters->q;
      pri_omega = 0.5*( current_machines[0]->parameters->omega + current_machines[1]->parameters->omega);
    }


  spikecountmean[0] += (pri_q/2.0)*pri_omega;
  spikecountssd[0] += ((pri_q/2.0)*pri_omega) * ((pri_q/2.0)*pri_omega) ;
  spikecounttotal[0] += 1;

  spikecountmean[1] += (pri_q/2.0)*pri_omega;
  spikecountssd[1] += ((pri_q/2.0)*pri_omega) * ((pri_q/2.0)*pri_omega) ;
  spikecounttotal[1] += 1;

  spikecountmean[2] += (pri_q)*pri_omega;
  spikecountssd[2] += ((pri_q)*pri_omega) * ((pri_q)*pri_omega) ;
  spikecounttotal[2] += 1;

  spikecountmean[3] += (pri_q)*pri_omega;
  spikecountssd[3] += ((pri_q)*pri_omega) * ((pri_q)*pri_omega) ;
  spikecounttotal[3] += 1;

  baselineonlymean[0] += pri_q/2.0;
  baselineonlyssd[0] += (pri_q/2.0)*(pri_q/2.0);
  baselineonlytotal[0] += 1;

  baselineonlymean[1] += pri_q/2.0;
  baselineonlyssd[1] += (pri_q/2.0)*(pri_q/2.0);
  baselineonlytotal[1] += 1;

  baselineonlymean[2] += pri_q;
  baselineonlyssd[2] += pri_q * pri_q;
  baselineonlytotal[2] += 1;

  baselineonlymean[3] += pri_q;
  baselineonlyssd[3] += pri_q * pri_q;
  baselineonlytotal[3] += 1;

}
*/
void YStateAndCountStats::Initialize(vector<Machine*>& current_machines, vector<Machine*>& proposed_machines)
{
  vector<Machine*> machine_s1_s2_s12;
  int i, j, T;
  double tmp_lagstate;
  double tmp_exsq;
  double lag1, lag2, lag12;
  int y;
  int use2;
  double pri_q, pri_omega;
  spikecountmean.clear();
  spikecountssd.clear();
  spikecounttotal.clear();
  baselineonlymean.clear();
  baselineonlyssd.clear();
  baselineonlytotal.clear();
  spikesums.clear();
  spikecounts.clear();
  
 if(current_machines.size() == 1)
    {
      pri_q = current_machines[0]->parameters->q;
      pri_omega = current_machines[0]->parameters->omega;
    }
  else
    {
      pri_q = current_machines[0]->parameters->q + current_machines[1]->parameters->q;
      pri_omega = 0.5*( current_machines[0]->parameters->omega + current_machines[1]->parameters->omega);
    }

  if(current_machines.size() == 1)
    {
      machine_s1_s2_s12.push_back(proposed_machines[0]);
      machine_s1_s2_s12.push_back(proposed_machines[1]);
      machine_s1_s2_s12.push_back(current_machines[0]);
      use2 = 0;
     
    }
  else
    {
      machine_s1_s2_s12.push_back(current_machines[0]);
      machine_s1_s2_s12.push_back(current_machines[1]);
      machine_s1_s2_s12.push_back(proposed_machines[0]);
      use2 = 1;
    }
  
  for(i = 0 ; i < 3; i++)
    {
      max_lag_s1_s2_s12.push_back(0.0);
    }

  T = machine_s1_s2_s12[0]->parameters->T;
  for(i = 0; i < T; i++)
    {
      for(j = 0; j < 3; j++)
	{
	  tmp_lagstate = machine_s1_s2_s12[j]->parameters->Lagstate(i);
	  if(tmp_lagstate > max_lag_s1_s2_s12[j])
	    max_lag_s1_s2_s12[j] = tmp_lagstate;
	}
    }
  for(j = 0 ; j < 3; j++)
    {
      if(max_lag_s1_s2_s12[j] > baseline_lag)
	max_lag_s1_s2_s12[j] = baseline_lag;
    }
  
  for(i = 0; i < 4; i++)
    {
      spikecountmean.push_back(0.0);
      spikecountssd.push_back(0.0);
      spikecounttotal.push_back(0.0);
      baselineonlymean.push_back(0.0);
      baselineonlyssd.push_back(0.0);
      baselineonlytotal.push_back(0.0);
      spikesums.push_back(0.0);
      spikecounts.push_back(0.0);
    }

  //order is so os sd ds ss
  spikesums.push_back(0.0);
  spikecounts.push_back(0.0);

  //setting some priors

  for(i = 0; i < T; i++)
    {
      lag1 = machine_s1_s2_s12[0]->parameters->Lagstate(i);
      lag2 = machine_s1_s2_s12[1]->parameters->Lagstate(i);
      lag12= machine_s1_s2_s12[2]->parameters->Lagstate(i);
      if(use2)
	{
	  y = machine_s1_s2_s12[0]->data->GetTotalCount(i) + machine_s1_s2_s12[1]->data->GetTotalCount(i);
	}
      else
	{
	  y = machine_s1_s2_s12[2]->data->GetTotalCount(i);
	}

      if(lag12 == 0)
	{
	  spikecountmean[3] += y;
	  spikecountssd[3] += y*y;
	  spikecounttotal[3] += 1;
	}
      if(lag12 >= max_lag_s1_s2_s12[2])
	{
	  baselineonlymean[3] += y;
	  baselineonlyssd[3] += y*y;
	  baselineonlytotal[3] += 1;
	}

      if( lag2 == -1)
	{
	  if(lag1 == 0)
	    {
	      spikecountmean[0] += y;
	      spikecountssd[0] += y*y;
	      spikecounttotal[0] += 1;

	      spikesums[0] += (double)(y);
	      spikecounts[0] += 1.0;
	    }
	  else
	    {
	      if(lag1 >= max_lag_s1_s2_s12[0])
		{
		  baselineonlymean[0] += y;
		  baselineonlyssd[0] += y*y;
		  baselineonlytotal[0] += 1;		  
		}
	    }
	}
      if( lag1 == -1)
	{
	  if(lag2 == 0)
	    {
	      spikecountmean[1] += y;
	      spikecountssd[1] += y*y;
	      spikecounttotal[1] += 1;

	      spikesums[1] += (double)(y);
	      spikecounts[1] += 1.0;
	    }
	  else
	    {
	      if(lag2 >= max_lag_s1_s2_s12[1])
		{
		  baselineonlymean[1] += y;
		  baselineonlyssd[1] += y*y;
		  baselineonlytotal[1] += 1;		  
		}
	    }
	}
      if( (lag1 == 0) && (lag2 == 0))
	{
	  spikecountmean[2] += y;
	  spikecountssd[2] += y*y;
	  spikecounttotal[2] += 1;
	  spikesums[4] += (double)(y);
	  spikecounts[4] += 1.0;
	}
      if( (lag1 >= max_lag_s1_s2_s12[0]) && (lag2 >= max_lag_s1_s2_s12[1]))
	{
	  baselineonlymean[2] += y;
	  baselineonlyssd[2] += y*y;
	  baselineonlytotal[2] += 1;
	}
      if( (lag1 >= max_lag_s1_s2_s12[0]) && (lag2 == 0))
	{
	  spikesums[3] += (double)(y);
	  spikecounts[3] += 1.0;
	}
      if( (lag1 == 0) && (lag2 >= max_lag_s1_s2_s12[1]))
	{
	  spikesums[2] += (double)(y);
	  spikesums[2] += 1.0;
	}

    }

  //now calculate the means
  for(i = 0; i < 4; i++)
    {
      if(spikecounttotal[i] > 0)
	{
	  if(spikecounttotal[i] > 1)
	    {
	      spikecountssd[i] = ((double)(spikecounttotal[i])*spikecountssd[i] - (spikecountmean[i]*spikecountmean[i]))/((double)(spikecounttotal[i])*((double)(spikecounttotal[i])-1.0));
	      spikecountssd[i] = sqrt(spikecountssd[i]);
	    }
	  else
	    {
	      spikecountssd[i] = 0.0;
	    }
	  spikecountmean[i] /= (double)(spikecounttotal[i]);
	}
      if(baselineonlytotal[i] > 0)
	{
	  if(baselineonlytotal[i] > 1)
	    {
	      baselineonlyssd[i] = ((double)(baselineonlytotal[i])*baselineonlyssd[i] - (baselineonlymean[i]*baselineonlymean[i]))/((double)(baselineonlytotal[i])*((double)(baselineonlytotal[i])-1.0));
	      baselineonlyssd[i] = sqrt(baselineonlyssd[i]);
	    }
	  else
	    {
	      baselineonlyssd[i] = 0.0;
	    }
	  baselineonlymean[i] /= (double)(baselineonlytotal[i]);
	}
      
    }
  
}


#include "../datparobjects.h"

/*=====================Data Objects (For Networks)==========================*/

/* Copy constructor */
netw_dat_obj::netw_dat_obj(netw_dat_obj* cpy)
 {
    
    int i;
    if(cpy->is_allocated)
      {
	/*fl = new int[cpy->T];*/
	/*ip = new int[cpy->T];*/
	for(i = 0; i < cpy->T; i++)
	  {
	    fl.push_back(cpy->fl[i]);
	    ip.push_back(cpy->ip[i]);
	  }
      }
    T = cpy->T;
    timezone = cpy->timezone;
    beginval = cpy->beginval;
    endval = cpy->endval;
    tstar = cpy->tstar;
    row_id = cpy->row_id;
    block_id = cpy->block_id;
    is_allocated = cpy->is_allocated;	       
    total_objects = cpy->total_objects;
  }

/* File object pointers for processing many lines of input in a row */
netw_dat_obj::netw_dat_obj(int N, int rid, ifstream &spikeinfo, ifstream &flinfo, ifstream &ipinfo)
  {
    /*RID is 0-indexed!*/
    /*For processing many lines in a row*/
    is_allocated = 0;
    GetData(N, rid, spikeinfo, flinfo, ipinfo);
  }

/*For obtaining data from just one file stream (uses head/tail command line) */
netw_dat_obj::netw_dat_obj(int N, int rid)
 {
    /*RID is 0-indexed!*/
    /*For obtaining data from just one file stream (uses head/tail command line) */
    int flipline;
    int infoline;

    is_allocated = 0;
    flipline = rid + 1;
    infoline = rid + 2; /*accounts for titles in file*/
    /* cout << "IP/FL line to read="<< flipline <<", info line to read=" << infoline << endl;
       cout.flush();*/
    string startcommand = "head -";
    string endipcommand = " " + ippath_name + " | tail -1 > _tmp_ip.txt";
    string endflcommand = " " + flpath_name + " | tail -1 > _tmp_fl.txt";
    string endinfocommand = " " + inpath_name + " | tail -1 > _tmp_info.txt";

    string ipcomm = startcommand + to_string(flipline) + endipcommand;
    string flcomm = startcommand + to_string(flipline) + endflcommand;
    system(ipcomm.c_str());
    system(flcomm.c_str());
         
    string infocomm = startcommand + to_string(infoline) + endinfocommand;
    system(infocomm.c_str());
   
    ifstream spinfo("_tmp_info.txt"); 
    ifstream ipinfo("_tmp_ip.txt");
    ifstream flinfo("_tmp_fl.txt");
    
    GetData(N, rid, spinfo, flinfo, ipinfo);

    ipinfo.close();
    flinfo.close();
    spinfo.close();
    
    if(!KEEP_TEST)
      {
	unlink("_tmp_ip.txt");
	unlink("_tmp_fl.txt");
	unlink("_tmp_info.txt");
      }
 }

netw_dat_obj::~netw_dat_obj()
  {
    /*
    if(is_allocated)
      {
	delete fl;
	delete ip;
      }
      is_allocated =0;
    */
  }

void netw_dat_obj::Print(ostream& str, int N)
{
 int i;
 set<machine_id_t>::iterator machine_iter;
 if((N == -1) || (N > T))
   {
     N = T;
   }
 str << "**" << endl;
 str << "Information for data object corresponding to row " << row_id << ": (" << T << " observations)" <<  endl;
 str << "Begin End TZ Tstar TotalObjects BlockID" << endl << beginval << " " << endval << " " << timezone << " " << tstar << " " << total_objects << " " << block_id << endl;
 str << "First " << N << " IP counts: "; 
 for(i = 0; i < (N-1); i++)
   str << " " << ip[i];
 str << endl;
 str << endl;
 str << "First " << N << " Flow counts: ";
 for(i = 0; i < N; i++)
   str << " " << fl[i];
 str << endl << "Associated Machine IDs:";
 for(machine_iter = machines.begin(); machine_iter != machines.end(); machine_iter++)
   str << (*machine_iter)  << " ";
 str << endl;
 str << endl<< "**" << endl;
 return;
}


void netw_dat_obj::GetData(int N, int rid, ifstream &spikeinfo, ifstream &flinfo, ifstream &ipinfo)
  {
    int i, tmp1, tmp2;
    /*fl = new int[N];
      ip = new int[N];*/
    is_allocated =1;
    T = N;   
    row_id = rid;
    spikeinfo >> timezone >> tstar >> beginval >> endval >> total_objects >> block_id;
    //cout << timezone << " " << tstar << " " << beginval << " " << endval << endl;
    for(i = 0; i < N; i++)
      {
	flinfo >> tmp1;
        fl.push_back(tmp1);
	ipinfo >> tmp2;
	ip.push_back(tmp2);
	
	/*cout << i << " " << fl[i] << " " << ip[i] << ":";*/

	}   
    return;
  }



/*=====================Data Objects (For Machines)==========================*/

/*Network blocks*/
NetworkBlock::NetworkBlock(NetworkBlock* cpy)
{
  block_allocations = cpy->block_allocations;
}
NetworkBlock::NetworkBlock(unsigned int uid, vector<short int>& indexes, vector<int>& totals)
{
  int i;
  int added = 0;
  if(indexes.size() == totals.size())
    {
      for(i = 0; i < totals.size(); i++)
	if(totals[i] > 0)
	  {
	    block_allocations[uid][indexes[i]] = totals[i];
	    added = 1;
	  }
    }						
  else
    cout << "Error: Network Block constructor: index and totals must be the same length" << endl;
  if(!added)
    block_allocations[uid] = network_inner_t();
}


NetworkBlock::NetworkBlock(unsigned int uid, vector<int>& totals)
{
  int i;
  int added =0;
  for(i = 0; i < totals.size(); i++)
    if(totals[i] > 0)
      {
	block_allocations[uid][i] = totals[i];
	added = 1;
      }
  if(!added)
    block_allocations[uid] = network_inner_t();
}



void NetworkBlock::Print(ostream& str, int N)
{
  network_outer_t::iterator outer_loop;
  network_inner_t::iterator inner_loop;
  for(outer_loop = block_allocations.begin(); outer_loop != block_allocations.end(); outer_loop++)
    {
      str << outer_loop->first << ":";
      for(inner_loop = outer_loop->second.begin(); inner_loop != outer_loop->second.end(); inner_loop++)
	{
	  if( (inner_loop->first < N) || (N < 0))
	    str << "(" << inner_loop->first << "," << inner_loop->second << ")";
	}
      str << endl;
    }
  return;
}

void NetworkBlock::Add(NetworkBlock* a)
{
  network_outer_t::iterator outer_loop;
  network_inner_t::iterator inner_loop;
  network_outer_t::iterator current_outer;
  network_inner_t::iterator current_inner;
  for(outer_loop = a->block_allocations.begin(); outer_loop != a->block_allocations.end(); outer_loop++)
    {
      current_outer = block_allocations.find(outer_loop->first);
      if(current_outer == block_allocations.end())
	{
	  //the network does not exist at all so just add the whole thing
	  block_allocations.insert(current_outer, *(outer_loop));
	}
      else
	{
	  //pairwise add elements to an existing network
	  for(inner_loop = outer_loop->second.begin(); inner_loop != outer_loop->second.end(); inner_loop++)
	    {
	      current_inner = current_outer->second.find(inner_loop->first);
	      if(current_inner == current_outer->second.end())
		{
		  //the time point does not exist; add it in
		  current_outer->second.insert(current_inner, *(inner_loop));
		}
	      else
		{
		  //add the value to the existing one
		  current_inner->second += inner_loop->second;
		}
	    }
	}
    }
  
}

int NetworkBlock::Subtract(NetworkBlock* a)
{
  //make sure that (this) is >= a at all points first
  network_outer_t::iterator outer_loop;
  network_inner_t::iterator inner_loop;
  network_outer_t::iterator current_outer;
  network_inner_t::iterator current_inner;
  network_outer_t tmp_allocations = block_allocations;
  int valid = 1;
  
  outer_loop = a->block_allocations.begin();
  while( (outer_loop != a->block_allocations.end()) && (valid))
    {
      current_outer = tmp_allocations.find(outer_loop->first);
      if(current_outer == tmp_allocations.end())
	{
	  //unless all entries for that one are zero it is invalid
	  if(!(a->isNonZeroNetwork(outer_loop->first)))
	     valid = 0;
	}
      else
	{
	  inner_loop = outer_loop->second.begin();
	  while( (inner_loop != outer_loop->second.end()) &&(valid))
	    {
	      current_inner = current_outer->second.find(inner_loop->first);
	      if( current_inner == current_outer->second.end())
		{
		  if(inner_loop->second > 0)
		    {
		      //subtracting a positive number from 0 in the current allocation is not allowed
		      valid = 0;
		    }
		}
	      else
		{
		  if(current_inner->second < inner_loop->second)
		    {
		      valid = 0;
		    }
		  else
		    {
		      current_inner->second -= inner_loop->second;
		      if(current_inner->second == 0)
			{
			  //if we've subtracted all of the flows, remove that entry from the inner loop
			  current_outer->second.erase(current_inner);
			}
		    }
		}
	      inner_loop++;
	    }
	  //if we've removed all of the flows from an entire network, set that network to the Zero network
	
	}
      outer_loop++;
    }
  if(valid)
    block_allocations = tmp_allocations;
  else
    {
      cout << "Error: Subtraction of block allocations would result in negative values: subtraction is not performed" << endl;
    }
  return valid;
}
 
void NetworkBlock::AddNetwork(unsigned int uid)
{
  if(block_allocations.find(uid)  == block_allocations.end())
    block_allocations[uid] = network_inner_t();
}

int NetworkBlock::GetCount(unsigned int uid, short int t)
{
  int count = 0;
  if(block_allocations.find(uid) != block_allocations.end())
    {
      if(block_allocations[uid].find((short int)t) != block_allocations[uid].end())
	count = block_allocations[uid][t];
    }
  return count;
}

void NetworkBlock::SetCountOnExistingNetwork(unsigned int uid, short int t, int val)
{
  if(block_allocations.find(uid) != block_allocations.end())
    {
      block_allocations[uid][t] = val;
    }
}

int NetworkBlock::RemoveZeroNetwork(unsigned int uid)
{
  int removed = 0;
  if(!(isNonZeroNetwork(uid)))
    {
      block_allocations.erase(uid);
      removed = 1;
    }
  return removed;
}

void NetworkBlock::ListNetworks(vector<unsigned int>* network_list)
{
  network_outer_t::iterator outer_loop;
  network_list->clear();
  for(outer_loop = block_allocations.begin(); outer_loop != block_allocations.end(); outer_loop++)
    {
      network_list->push_back(outer_loop->first);
    }
}

int NetworkBlock::NumNetworks()
{
  return block_allocations.size();
}

int NetworkBlock::isNonZeroNetwork(int netid)
  {
    network_outer_t::iterator outer_network;
    network_inner_t::iterator inner_loop;
    int nonzero_entry_exists = 0;
    outer_network = block_allocations.find(netid);
    if(outer_network != block_allocations.end())
      {
	inner_loop = outer_network->second.begin();
	while( (inner_loop != outer_network->second.end()) && (!nonzero_entry_exists))
	  {
	    if(inner_loop->second > 0)
	      nonzero_entry_exists = 1;
	    inner_loop++;
	  }
      }
    return nonzero_entry_exists;
  }

void NetworkBlock::SumBlocks(vector<int>& totals)
{
  int i;
  network_outer_t::iterator outer_loop;
  network_inner_t::iterator inner_loop;
  for(i = 0; i < totals.size(); i++)
    totals[i] = 0;
  for(outer_loop = block_allocations.begin(); outer_loop != block_allocations.end(); outer_loop++)
    {
      for(inner_loop = outer_loop->second.begin(); inner_loop != outer_loop->second.end(); inner_loop++)
	{
	  if(inner_loop->first < totals.size())
	    totals[inner_loop->first] += inner_loop->second;
	  else
	    cout << "Error in NetworkBlock::SumBlocks(vector<int>& totals) : Block allocations contain index greater than the length of totals vector" << endl;
	}
    }
}

/* Copy constructor */
mach_dat_obj::mach_dat_obj(mach_dat_obj* cpy)
 {
    int i;
    if(cpy->is_allocated)
      {
	for(i = 0; i < cpy->T; i++)
	  {
	    fl.push_back(cpy->fl[i]);
	  }
      }
    T = cpy->T;
    beginval = cpy->beginval;
    endval = cpy->endval;
    timezone = cpy->timezone;
    tstar = cpy->tstar;
    is_allocated = cpy->is_allocated;
    networks = new NetworkBlock(cpy->networks);
    fl_nb_agree = cpy->fl_nb_agree;
    uid = cpy->uid;
    iteration_born = cpy->iteration_born;
  }

/*nearly copy constructor*/
mach_dat_obj::mach_dat_obj(mach_dat_obj* cpy, machine_id_t id, int iteration)
 {    
    int i;
    if(cpy->is_allocated)
      {
	/*fl = new int[cpy->T];*/
	/*ip = new int[cpy->T];*/
	for(i = 0; i < cpy->T; i++)
	  {
	    fl.push_back(cpy->fl[i]);
	  }
      }
    T = cpy->T;
    beginval = cpy->beginval;
    endval = cpy->endval;
    timezone = cpy->timezone;
    tstar = cpy->tstar;
    is_allocated = cpy->is_allocated;
    networks = new NetworkBlock(cpy->networks);
    fl_nb_agree = cpy->fl_nb_agree;
    uid = id;
    iteration_born = iteration;
  }


mach_dat_obj::mach_dat_obj(netw_dat_obj* cpy, machine_id_t id, int iteration)
 {
    
    int i;
    if(cpy->is_allocated)
      {
	/*fl = new int[cpy->T];*/
	/*ip = new int[cpy->T];*/
	for(i = 0; i < cpy->T; i++)
	  {
	    fl.push_back(cpy->fl[i]);
	  }
      }
    T = cpy->T;
    beginval = cpy->beginval;
    endval = cpy->endval;
    timezone = cpy->timezone;
    tstar = cpy->tstar;
    is_allocated = cpy->is_allocated;
    networks = new NetworkBlock(cpy->block_id, fl); 
    fl_nb_agree = 1;
    uid = id;
    iteration_born = iteration;
  }

void mach_dat_obj::SumBlocks()
{
  if(fl_nb_agree == 0)
    {
      networks->SumBlocks(fl);
      fl_nb_agree = 1;
    }
}

int mach_dat_obj::NumNetworks()
{
  return(networks->NumNetworks());
}

int mach_dat_obj::GetTotalCount(int t)
{
  int count =0;
  SumBlocks();
  if( (t >= 0) && (t < fl.size()))
    count = fl[t];
  else
    cout << "Error, requesting count at time " << t << " greater than data size " << fl.size() << endl;
  return count;
}

int mach_dat_obj::GetSubCount(unsigned int networkid, int t)
{
  return networks->GetCount(networkid, t);
}

void mach_dat_obj::SetCountExisting(unsigned int networkid, int hour, int count)
{
  networks->SetCountOnExistingNetwork(networkid, hour, count);
  fl_nb_agree = 0;
}

int mach_dat_obj::SubtractBlocks(mach_dat_obj* machine)
{
  int success;
  success = networks->Subtract(machine->networks);
  if(success)
    fl_nb_agree = 0;
  return(success);
}

void mach_dat_obj::AddBlocks(mach_dat_obj* machine)
{
  networks->Add(machine->networks);
  fl_nb_agree = 0;
}

void mach_dat_obj::ListNetworks(vector<unsigned int>* network_list)
{
  networks->ListNetworks(network_list);
}


void mach_dat_obj::AddNetwork(unsigned int networkid)
{
  networks->AddNetwork(networkid);
}

int mach_dat_obj::RemoveZeroNetwork(unsigned int networkid)
{
  int success;
  success = networks->RemoveZeroNetwork(networkid);
  return(success);
}

mach_dat_obj::~mach_dat_obj()
  {
    delete networks;
  }

void mach_dat_obj::Print(ostream& str, int N)
{
 int i;
    if((N == -1) || (N > T))
      {
	N = T;
      }
    str << "*****" << endl;
    str << "Begin End TZ Tstar" << endl << beginval << " " << endval << " " << timezone << " " << tstar << " " << endl;
    str << endl;
    str << endl;
    str << "First " << N << " Flow counts: ";
    for(i = 0; i < N; i++)
      str << " " << fl[i];
    str << endl << "Block allocations: " << endl;
    networks->Print(str,N);
    str << endl<< "*****" << endl;
    return;
}

/*=====================Parameter Object (For a single Machine)==========================*/

curr_par_obj::curr_par_obj()
  {
    is_allocated = 0;
    T = 0;
    tstar = 0;
    birthdate = 0;
    deathdate = 0;
    is_immune = 0;
    num_off_off = 0; 
    num_off_state_seq=0;
 }


/* Set current parameters with starting state for poisson and transition, and likely state values given observed data structure in dat_obj object data */
curr_par_obj::curr_par_obj(mach_dat_obj* data, double* psi_params, double q_start, double alpha_start, double omega_start, double* rho_start, double* gamma_start, double * nu_start, double off_min, double off_max, int immune)
  {
    int i;
    T = data->T;
    for(i = 0; i < T; i++)
      {
	state.push_back(0);
	lagstate.push_back(0.0);
      }
    /*    state = new short int[T];
	  lagstate = new double[T];*/
    is_allocated = 1;
    tstar = data->tstar;
    birthdate = 0;
    deathdate = T;
    is_immune = immune;
    num_off_off = 0;
    num_off_state_seq = 0;
    q = q_start;
    alpha = alpha_start;
    omega = omega_start;

    for(i = 0; i<=2; i++)
      {
	rho[i] = rho_start[i];
	gamma[i] = gamma_start[i];
	nu[i] = nu_start[i];
      }

    SetPsiParameterization(psi_params);
    SetStateStart(data, off_min, off_max);
    SetLagstate();  
  }

/*"sort of" copy constructor (copy parameters but initialize state)*/
curr_par_obj::curr_par_obj(mach_dat_obj* data, curr_par_obj* cpy, hype_par_obj* hyper)
{
  int i,j,k;
    T = data->T;
    for(i = 0; i < T; i++)
      {
	state.push_back(0);
	lagstate.push_back(0.0);
      }
    /*    state = new short int[T];
	  lagstate = new double[T];*/
    is_allocated = 1;
    tstar = data->tstar;
    birthdate = 0;
    deathdate = T;
    is_immune = cpy->is_immune;
    num_off_off = 0;
    num_off_state_seq = 0;
    q = cpy->q;
    alpha = cpy->alpha;
    omega = cpy->omega;

    for(i = 0; i<=2; i++)
      {
	rho[i] = cpy->rho[i];
	gamma[i] = cpy->gamma[i];
	nu[i] = cpy->nu[i];
      }

    for(i = 0; i < 3; i++)
      for(j = 0; j < 3; j++)
	for(k = 0; k < 2; k++)
	  psi_parameterization[i][j][k] = cpy->psi_parameterization[i][j][k];
    SetStateStart(data, hyper->off_min, hyper->off_max);
    SetLagstate();  
}

curr_par_obj::curr_par_obj(ifstream& init_pars)
{
  int i, j,k;
  is_allocated = 0;
  T = 0;
  tstar = 0;
  for(i = 0; i < 3; i++)
    for(j =0; j < 3; j++)
      for(k=0; k < 2; k++)
	init_pars >> psi_parameterization[i][j][k];
  init_pars >> q >> alpha >> omega;
  for(i =0; i < 3; i++)
    init_pars >> rho[i];
  for(i = 0; i < 3; i++)
    init_pars >> gamma[i];
  for(i = 0; i < 3; i++)
    init_pars >> nu[i];
  is_immune = 0;
  birthdate = 0;
  deathdate = 0;
  num_off_off = 0;
  num_off_state_seq = 0;
}


/*Copy constructor */
curr_par_obj::curr_par_obj(curr_par_obj* cpy)
  {
    int i;
    if(cpy->is_allocated)
      {
	T = cpy->T;
	for(i = 0 ; i< T; i++)
	  {
	    state.push_back(0);
	    lagstate.push_back(0.0);
	  }
	is_allocated = 1;
	InternalSet(cpy);
      }
    else
      {
         is_allocated = 0;
	 T = 0;
	 is_immune = 0;
	 num_off_off = 0;
	 birthdate = 0;
	 deathdate = 0;
	 num_off_state_seq = 0;
	 tstar = 0;
      }
  }


void curr_par_obj::SetFromCopy(curr_par_obj* cpy)
 {
   T=cpy->T; 
   while(lagstate.size() > T)
     lagstate.pop_back();
   while(state.size() > T)
     state.pop_back();
   
   while(lagstate.size() < T)
     lagstate.push_back(0.0);
   while(state.size() < T)
     state.push_back(0);

    is_allocated = 1;
    InternalSet(cpy);
  }

void curr_par_obj::SetState(int index, short int val)
  {
    double newlag;
    int i;
    int oldstate;
    oldstate = state[index];
    if(val != oldstate)
      {	
	state[index] = val;
	AdjustSingleOffState(index, oldstate, val);
	i = index;
	newlag = SetIndividualLagState(i);
	while( (i < T) && (newlag != lagstate[i]))
	  {
	    lagstate[i] = newlag;
	    i++;
	    newlag = SetIndividualLagState(i);
	  }
      }
  }

void curr_par_obj::AdjustSingleOffState(int index, short int oldstate, short int newstate)
{
  int addsub = 0;
  int countval = 0;
  int check_previous;
  int check_following;
  int state_trans_add = 0;
  if( (index <= deathdate) && (index >= birthdate))
    {
      check_previous = !(index == birthdate);
      check_following = !(index == deathdate);
      if(!(oldstate == newstate))
	{
	  if(check_previous)
	    {
	      countval += (state[index-1] == 0);
	    }
	  if(check_following)
	    {
	      countval += (state[index+1] == 0);
	    }
	  if(oldstate == 0)
	    {
	      addsub = -1;
	      //newstate is 1 or 2
	      //if 2 0's, breaking up 1 run into 2 (+1)
	      //if 1 0, no change (0)
	      //if 0 0's, subtracting a run of only one 0 (-1)
	    }
	  if(newstate == 0)
	    {
	      addsub = 1;
	      //oldstate is 1 or 2
	      //if 2 zeroes, making 2 runs into one (-1)
	      //if 1 zero, no change (0)
	      //if 0 zeroes adding a run of only one (+1)
	    }
	}
      num_off_off += (addsub)*(countval);
      num_off_state_seq += (-1*addsub)*(countval-1);
      if( (num_off_off > (deathdate - birthdate)) || (num_off_off < 0))
	{
	  cout << "Error in setting number of off->off transitions, number =" << num_off_off << endl;
	}
    }
  return;
}

void curr_par_obj::SetBirthDeathImmune(short int birth, short int death, short int immune)
{
  short int oldbirth = birthdate;
  short int olddeath = deathdate;
  short int oldimmune = is_immune;
  short int old_end_trunc;
  short int new_end_trunc;
  int is_valid =1;
  int adj = 0;
  int birthsub = 0;
  int num_birth_transitions = 0;
  
  int adj2 = 0;
  int deathsub = 0;
  int num_death_transitions = 0;
  if(is_allocated)
    {
      if(olddeath >= T)
	old_end_trunc = T-1;
      else
	old_end_trunc = olddeath;
      
      if(death >= T)
	new_end_trunc = T-1;
      else
	new_end_trunc = death;
   
      if( (immune == 1) && (death < T))
	{
	  is_valid = 0;
	}
      if(birth != oldbirth)
	{
	  if(birth > oldbirth)
	    {
	      adj = (state[birth] != 0);
	      num_birth_transitions = birth - oldbirth - adj;
	      birthsub = -1;
	    }
	  else
	    {
	      adj = (state[oldbirth] != 0);
	      num_birth_transitions = oldbirth - birth - adj;
	      birthsub = 1;
	    }
	}
      
      if(death != olddeath)
	{
	  if(death < olddeath)
	    {
	      adj2 = (state[new_end_trunc] != 0);
	      num_death_transitions = old_end_trunc - new_end_trunc - adj2;
	      deathsub = -1;
	    }
	  else
	    {
	      adj2 = (state[old_end_trunc] != 0);
	      num_death_transitions = new_end_trunc - old_end_trunc - adj2;
	      deathsub = 1;
	    }
	}
      
      if(is_valid)
	{
	  is_immune = immune;
	  birthdate = birth;
	  deathdate = death;
	  num_off_off += (birthsub*num_birth_transitions) + (deathsub*num_death_transitions);
	  num_off_state_seq += (birthsub*adj) + (deathsub*adj2);
	  if( (num_off_off > (deathdate - birthdate)) || (num_off_off < 0))
	    {
	      cout << "Error in setting number of off->off transitions, number =" << num_off_off << endl;
	    }
	}
      else
	{
	  cout << "Error in SetBirthDeathImmune: Cannot set immune = 1 and Deathdate = " << death << " < Max time" << endl;
	}
    }
  else
    {
      birthdate = birth;
      deathdate = death;
      is_immune = immune;
    }
  return;
}

int curr_par_obj::BDIindex()
{
  int subind = 0;
  int index;
  index = 4*(birthdate != 0) + 2*(deathdate == T) + (is_immune == 1);
  if(index > 0)
    subind++;
  if(index > 5)
    subind++;
  index = index - subind;
  return(index);
}


void curr_par_obj::SetPsiParameterization(double* pars)
  {
    int i,j,k,z;
    z = 0;
    for(i = 0; i <=2; i++)
      for(j = 0; j <=2; j++)
	for(k = 0; k <=1; k++)
	  {
	    psi_parameterization[i][j][k] = pars[z];
	    z++;
	  }
  }  



 void curr_par_obj::SetStateStart(mach_dat_obj* data, double off_min, double off_max)
  {
    
    int nonzero_count;
    int i,j, end_trunc;
    double topNcutoff;
    short int* topN;
    int ind;
    double delta;
    int last_spike, decay_zero_count;
    int in_off_seq;
    int num_off_check = 0;
    num_off_off = 0;
    num_off_state_seq = 0;
    nonzero_count = 0;
    in_off_seq = 0;
    birthdate = 0;
    deathdate = T;
    end_trunc = deathdate - (deathdate == T);
    for(i = 0; i < T; i++)
      {
	if(data->GetTotalCount(i) > 0)
	  nonzero_count++;	
      }
    if(nonzero_count > 1)
      {

	ind = (int) floor( (double)(nonzero_count - 1)*0.9);
	delta = (double)(nonzero_count -1 )*0.9 - (double)ind;
	topN = new short int[nonzero_count];

	j = 0;
	for(i =0; i < T; i++)
	  {
	    if(data->GetTotalCount(i) > 0)
	      {
		topN[j] = data->GetTotalCount(i);
		j++;
	      }
	  }
	/*get the 95th percentile of the non-zero data for the spike cutoff */
	gsl_heapsort(topN, nonzero_count, sizeof(short int), compare_shortint );
	topNcutoff = (1 - delta)*(double)topN[ind] + delta*(double)topN[ind + 1];
	/*	cout << "Top N Cutoff = " << topNcutoff << endl;
	cout << "topN = " << endl;
	for(i = 0; i < nonzero_count; i++)
	  cout << topN[i] << " ";
	  cout << endl;*/
	delete topN;
	last_spike = 0;
	decay_zero_count = 3;
	for(i = 0; i < T; i++)
	  {
	    if(data->GetTotalCount(i) == 0)
	      {
		if(decay_zero_count > 2)
		  {
		    state[i] = 0;
		    if( (i >= birthdate) && (i <= end_trunc))
		      {
			num_off_off++;
			if(!in_off_seq)
			  in_off_seq = 1;
		      }
		  }
		else
		  {
		    state[i] = 2;
		    decay_zero_count++;
		    if( (i >= birthdate) && (i <= end_trunc))
		      {
			if(in_off_seq)
			  {
			    in_off_seq = 0;
			    num_off_state_seq++;
			    num_off_off--;
			  }
		      }
		  }
		last_spike = 0;
	      }
	    else
	      {
		if( (i >= birthdate) && (i <= end_trunc))
		  {
		    if(in_off_seq)
		      {
			in_off_seq = 0;
			num_off_state_seq++;
			num_off_off--;
		      }
		  }
		decay_zero_count = 0;
		if(((double)(data->GetTotalCount(i)) > topNcutoff) && (!last_spike))
		  {
		    state[i] = 1;
		    last_spike = 1;
		  }
		else
		  {
		    state[i] = 2;
		    last_spike = 0;
		  }

	      }
	  }

      }
    else
      {
	if(nonzero_count == 0)
	  for(i = 0; i < T; i++)
	    {
	      if(data->GetTotalCount(i) == 0)
		{
		  state[i] = 0;
		  num_off_off++;
		  if(!in_off_seq)
		    in_off_seq =1;
		}
	      else
		{
		  state[i] = 2;
		  if(in_off_seq)
		    {
		      num_off_state_seq++;
		      in_off_seq = 0;
		      num_off_off--;
		    }
		}
	    }
      }
    if(in_off_seq)
      {
	in_off_seq =0;
	num_off_state_seq++;
	num_off_off--;
      }
    //cout << "num_off_off=" << num_off_off << " num_off_seq_count=" << num_off_seq_count << endl; 
    if(num_off_state_seq > 0)
      off_lambda = (float)(num_off_off)/(float)(num_off_state_seq);
    else
      off_lambda = off_min;
    if(off_lambda > off_max)
      off_lambda = off_max;
    if(off_lambda < off_min)
      off_lambda = off_min;

    
    return;
  }

double curr_par_obj::SetIndividualLagState(int i)
  {
    double current_lag, newlag;
    if(i == 0)
      current_lag = INFINITY;
    else
      current_lag = lagstate[i-1];
    
    if(state[i] == 0) //off
      {
	newlag = -1.0;     
      }
    else if(state[i] == 1) //spike
      {
	newlag = 0.0;
      }
    else if(state[i] == 2)
      {
	if( (current_lag == INFINITY) || (current_lag == -1.0))
	  newlag = INFINITY;
	else
	  newlag = current_lag + 1.0;
      }
    return(newlag);
  }
 
void curr_par_obj::SetLagstate()
  {
    int i;
    for( i = 0; i < T; i++)
      {
	lagstate[i] = SetIndividualLagState(i);
      }
  }


void curr_par_obj::Print(ostream& str, int N)
  {
    if( (N == -1) || (N > T))
      {
	N = T;
      }
    int i,j,k;
    str << "Parameter Values" << endl;
    str << "q="<< q << " alpha=" << alpha << " omega=" << omega << " offlambda=" << off_lambda;
    for(i = 0; i < 3; i++)
      {
	str << endl << "state " << i << ": rho=";
	str << rho[i] << " gamma=" << gamma[i] << " nu=" << nu[i] << " ";
      }
    str << endl;
    str << "Immune = " << is_immune << endl;
    str << "state vector ( birth = " << birthdate << ", death = " << deathdate << " ): "<<endl;
    for(i = 0; i < N; i++)
      str << state[i] << " ";
    str << endl << "lagstate:" << endl;
    for(i = 0; i < N; i++)
      str << lagstate[i] << " ";
    str << endl;
    str << "Psi parameterization: " << endl;
    for(i = 0; i <=2; i++)
      {
	str << i << ":" << endl;
	for(j = 0; j <=2; j++)
	  {
	    str << "  " << j << ":"; 
	    for(k = 0; k <=1; k++)
	      {
		str << psi_parameterization[i][j][k] << " ";
	      }
	    str << endl;
	  }
      }
    return;
  }


void curr_par_obj::PrintState(ostream &fileextension)
  {
    int i;
    for(i = 0; i < T; i++)
      {
	fileextension << state[i] << " ";
      }
    fileextension << endl;

  }


void curr_par_obj::PrintPars(ostream &fileextension, vector<int>* is_estimated)
  {
    int i;
    double parlist[16];
    parlist[0]=off_lambda;
    parlist[1]=q;
    parlist[2]=alpha;
    parlist[3]=omega;
    parlist[4] = rho[0];
    parlist[5] = gamma[0];
    parlist[6] = nu[0];
    parlist[7] = rho[1];
    parlist[8] = gamma[1];
    parlist[9] = nu[1];
    parlist[10]= rho[2];
    parlist[11]= gamma[2];
    parlist[12]= nu[2];

    parlist[13]= birthdate;
    parlist[14]= deathdate;
    parlist[15]= is_immune;
  
    if(is_estimated == NULL)
      {
	for(i = 0; i < 15; i++)
	  {
	    fileextension << parlist[i] << " ";
	  }
      }
    else
      {
	for(i = 0; i < 13; i++)
	  {
	    if(is_estimated->operator[](i) == 1)
	      fileextension << parlist[i] << " ";
	  }
	if(is_estimated->operator[](14) == 1)
	  {
	    fileextension << parlist[13] << " " << parlist[14] << " " << parlist[15];
	  }
      }
    fileextension << endl;
  }
    
double curr_par_obj::GetRate(int t)
{
  double rate = 0.0;
  if( (t >= birthdate) && ( t <= deathdate))
    {
      if(lagstate[t] == INFINITY)
	{
	  rate = q;
	}
      else
	if(lagstate[t] == -1)
	  {
	    rate = 0.0;
	  }
	else
	  {
	    rate = q*(1 + pow(alpha,lagstate[t])*(omega - 1.0));
	  }

    }
  return rate;
}


curr_par_obj::~curr_par_obj()
  {
    /*
    if(is_allocated)
      {
	delete state;
	delete lagstate;
      }
    */
 }

void curr_par_obj::InternalSet(curr_par_obj* cpy)
  {
    int i,j,k;
    tstar = cpy->tstar;
    birthdate = cpy->Birthdate();
    deathdate = cpy->Deathdate();
    is_immune = cpy->Is_immune();
    num_off_off = cpy->NumOffOffTrans();
    num_off_state_seq = cpy->NumOffSeq();
    T = cpy->T;
    for(i = 0; i < T; i++)
      {
	state[i] = cpy->State(i);
	lagstate[i] = cpy->Lagstate(i);
      }
    q = cpy->q;
    alpha = cpy->alpha;
    omega = cpy->omega;
    off_lambda = cpy->off_lambda;
    for(i = 0; i <= 2; i++)
      {
       	rho[i] = cpy->rho[i];
	gamma[i] = cpy->gamma[i];
	nu[i] = cpy->nu[i];
      }
    for(i = 0; i <=2; i++)
      for(j = 0; j <=2; j++)
	for(k = 0; k <=1; k++)
	  psi_parameterization[i][j][k] = cpy->psi_parameterization[i][j][k];
  }

int netw_par_obj::id = 0;

int netw_par_obj::getBoundaryID()
{
  return id++;
}

 

netw_par_obj::netw_par_obj()
{

  churn = 0.0;
  type = 0; /*singleton (0), NAT (1), DHCP (2) */
  boundary_id = getBoundaryID();
}

netw_par_obj::netw_par_obj(ifstream &netw_pars)
{
  netw_pars >> churn >> type;
  boundary_id = getBoundaryID();
}

void netw_par_obj::Print(ostream& str)
{
  str << "Network Parameters: churn="<<churn <<" type=" << type << " boundary_id="<< boundary_id << endl;
}

hype_par_obj::hype_par_obj(ifstream &hyper_par_file)
{
  hyper_par_file >> off_min >> off_max;
}

netw_par_obj::netw_par_obj(netw_par_obj* cpy)
{
  churn = cpy->churn;
  type = cpy->type;
  boundary_id = getBoundaryID();
}


global_par_obj::global_par_obj(ifstream &inputfile)
{
  int i;
  inputfile >> immune_rate >> immune_hyper[0] >> immune_hyper[1];
  for(i = 0; i < 5; i++)
    inputfile >> survival_hyper_mean_position[i];
  for(i = 0; i < 5; i++)
    inputfile >> survival_hyper_mean_height[i];
  inputfile >> survival_hyper_temperature;
  T = survival_hyper_mean_position[4]+1;
  survival_rate = new double[T];
  num_nonoff_nonimmune_transitions = new int[T];
  num_deaths = new int[T];
  //cout << "Survival Prior means" << endl;
  for(i = 0; i < T; i++)
    {
      num_nonoff_nonimmune_transitions[i] = 0;
      num_deaths[i] = 0;
      survival_rate[i] = SurvivalRatePriorMean(i);
      //cout << "(" << i << ","<< survival_rate[i] << ") ";
    }
  cout << endl;
  num_immune = 0;
  num_not_immune = 0;

}

double global_par_obj::SurvivalRatePriorMean(int t)
{
  //Return a piecewise linear prior dependent on the positions and heights
  int pos;
  int found, nonend;
  double meanval, slope;
  pos = 0;
  nonend = 0;
  found = 0;
  while(!found && (pos < 5))
    {
      if(t == survival_hyper_mean_position[pos])
	found = 1;
      else
	{
	  if(t < survival_hyper_mean_position[pos])
	    {
	      found = 1;
	      nonend = 1;
	    }
	  else
	    {
	      pos++;
	    }
	}
      
    }
  if(pos == 5)
    {
      cout << "Error in SurvivalRatePrior: time t=" << t << " out of range (0,"<<T-1<<")"<<endl;
    }
  else
    {
      if(nonend == 0)
	{
	  meanval = survival_hyper_mean_height[pos];
	}
      else
	{
	  slope = (survival_hyper_mean_height[pos] - survival_hyper_mean_height[pos-1])/( (double)survival_hyper_mean_position[pos] - (double)survival_hyper_mean_position[pos-1]);
	  meanval = survival_hyper_mean_height[pos-1] + (double)(t - survival_hyper_mean_position[pos-1])*slope;
	}
    }
  return meanval;
}

void global_par_obj::Print(ostream &str, int iteration_id)
{
  int i;
  if(iteration_id >= 0)
    str << iteration_id << " ";
  str << immune_rate << " ";
  for(i = 0; i < 5; i++)
    {
      str << survival_hyper_mean_position[i] << " " << survival_hyper_mean_height[i] << " ";
    }
  for(i = 0; i < T; i++)
    {
      str << survival_rate[i] << " ";
    }
  str << endl;
}

void global_par_obj::UpdateGlobalImmunityAndSurvivalCounts(curr_par_obj* parameters, int add_or_delete)
{
  int i;
  int nonoff_count_failure= 0;
  int end_trunc;
  if(parameters->Deathdate() >= parameters->T)
    end_trunc = parameters->T;
  else
    end_trunc = parameters->Deathdate();
  if(parameters->Is_immune())
    num_immune += add_or_delete;
  else
    {
      num_not_immune += add_or_delete;
      if(parameters->Deathdate() < parameters->T)
	num_deaths[parameters->Deathdate()] += add_or_delete;
      
      for(i = parameters->Birthdate(); i < end_trunc; i++)
	{
	  if(parameters->State(i) != 0)
	    {
	      num_nonoff_nonimmune_transitions[i] += add_or_delete;
	      if(num_nonoff_nonimmune_transitions[i] < 0)
		nonoff_count_failure = 1;
	    } 
	}
    }
  if( (num_deaths < 0) || (num_immune < 0) || (num_not_immune < 0) || (nonoff_count_failure < 0))
    cout << "Error in updating global parameters: counts < 0"<< endl;
}


global_par_obj::~global_par_obj()
{
  map<int, int*>::iterator bd;
  delete survival_rate;
  delete num_nonoff_nonimmune_transitions;
  delete num_deaths;
  for(bd = boundary_info.begin(); bd != boundary_info.end(); bd++)
    delete bd->second;
  boundary_info.clear();
}

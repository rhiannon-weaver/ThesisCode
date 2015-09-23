#include "../graphlinks.h"
#define ALPHA_SET 0.4
//index lower upper iterations usebirthdeath fix_decay initfile
int main(int argc, char* argv[])
{
  int lb = atoi(argv[2]);
  int ub = atoi(argv[3]);
  int iterations = atoi(argv[4]);
  string netdirectorypath = "/Users/rweaver/Research/analysis/Conficker/spikestatemodel/code/parameterfiles/datasets/SingleNets/Network";
  string initfilepath = "/Users/rweaver/Research/analysis/Conficker/spikestatemodel/code/parameterfiles/initfiles/init";
  string networkindex = argv[1];
  string suffix = ".txt";
  int usebirthdeath = atoi(argv[5]);
  int fix_decay = atoi(argv[6]);
  string initfile = argv[7];
  int i;
  
  string s = "/Users/rweaver/Research/analysis/Conficker/spikestatemodel/code/parameterfiles/hyperpar_initial_strong";
  string s2 = "_nobirthdeath";
  
  string par_init;
  if(usebirthdeath == 1)
     par_init = s+suffix;
  else
     par_init = s+s2+suffix;
  
  system(("python /Users/rweaver/Research/analysis/Conficker/spikestatemodel/code/testing/makeNetworkinitfiles.py " + networkindex + " " + initfile).c_str());

  NetworkGraph testgraph(1224, (netdirectorypath + networkindex).c_str(), "/Users/rweaver/Research/analysis/Conficker/spikestatemodel/code/parameterfiles/network_initial_pars.txt", "/Users/rweaver/Research/analysis/Conficker/spikestatemodel/code/parameterfiles/machine_initial_parameters_strong3.txt", par_init.c_str(),"parameterfiles/global_parameters.txt",lb, 0);
  testgraph.SetAllBoundaries(1);
  testgraph.InitializeMCMC((initfilepath + networkindex + suffix).c_str());
  if(usebirthdeath == 0)
    {
      testgraph.SetAllBirthDeathImmune(0,1224,1);
    }
  if(fix_decay == 1)
    {
      testgraph.SetAllDecayRates(ALPHA_SET);
    }
  testgraph.RunMCMC(iterations);

  for(i = lb+1; i <= ub; i++)
    {
      testgraph.ReinitializeMCMC(i,(initfilepath + networkindex + suffix).c_str(), 0);
      if(usebirthdeath == 0)
        {
	  testgraph.SetAllBirthDeathImmune(0,1224,1);
        }
      if(fix_decay)
	{
	  testgraph.SetAllDecayRates();
	}
      testgraph.RunMCMC(iterations);
    }  
  return 0;
}

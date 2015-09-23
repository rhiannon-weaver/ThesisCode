#include "../graphlinks.h"

int main(int argc, char* argv[])
{
int iter = atoi(argv[1]);
NetworkGraph testgraph(1224, "parameterfiles/datasets/Machine54", "parameterfiles/network_initial_pars.txt", "parameterfiles/machine_initial_parameters_strong3.txt", "parameterfiles/hyperpar_initial_strong.txt","parameterfiles/global_parameters.txt",1, 0);
testgraph.InitializeMCMC("parameterfiles/Machine54NonInformativeInit.txt");
testgraph.RunMCMC(iter);
return 0;
}

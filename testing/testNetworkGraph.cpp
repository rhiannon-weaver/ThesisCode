#include "../graphlinks.h"

int main()
{

Machine* selected_machine;
Machine* selected_machine_2;
Machine* merged_machine;
Network* selected_network;
vector<int> new_networks;
int i,j;

 NetworkGraph testgraph(1224, "parameterfiles/datasets/WaledacInformedSingletons", "parameterfiles/network_initial_pars.txt", "parameterfiles/machine_initial_parameters_strong3.txt", "parameterfiles/hyperpar_initial_strong.txt","parameterfiles/global_parameters.txt",0, 0);
testgraph.InitializeMCMC("parameterfiles/WaledacInitFile.txt");
testgraph.RunMCMC(50);

/*
for(i = 0; i < 20; i++)
{
  cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
//  testgraph.PrintNetworks(1);
  selected_network = testgraph.RandomlySelectNetwork();
//  selected_network->Print(20);

  selected_machine = testgraph.RandomlySelectMachine();
//  selected_machine->Print(20);
  testgraph.InsertMachine(new Machine(selected_machine));
  testgraph.DeleteMachine(selected_machine->uid);
}

for(i = 0; i < 10; i++)
{
  selected_machine = testgraph.RandomlySelectMachine();
  selected_machine_2 = testgraph.RandomlySelectMachine();

  merged_machine = new Machine(selected_machine);
  merged_machine->data->networks->Add(selected_machine_2->data->networks, new_networks);
  for(j = 0; j < new_networks.size(); j++)
	cout << new_networks[j] << " " ;
  cout << endl;
  while(new_networks.size() > 0)
    new_networks.pop_back();

  testgraph.DeleteMachine(selected_machine->uid);
  testgraph.DeleteMachine(selected_machine_2->uid);
  testgraph.InsertMachine(merged_machine);

  testgraph.PrintNetworks(1);
}
*/

return 0;

}

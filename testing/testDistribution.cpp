#include "../distributionlist.h"

int main()
{
  int i, k;
  DistributionList dists(0);
  double j;
 
  for(j = 0.0; j < 5000; j+=1.0) 
    cout << j <<" " <<  dists.LogDensity(4, j, 1.5, 3.0) << endl;

/*
  for(i = 0; i < 1000; i++)
	cout << j[i] << " ";
  cout << endl;
  dists.Simulate(4, 1000, j, 1.0, 2.0);
  for(i = 0; i < 1000; i++)
	cout << j[i] << " ";
  cout << endl;


  double pk[10];
  for(i = 0; i < 10; i++)
    { 
	pk[i] = (double)(i) /45.0;
	cout << pk[i] << " ";
    }
  cout << endl;

  unsigned int breakdown[10];
  for(i = 0; i < 100; i++)
  {
     dists.Simulate(8, 1, breakdown, 100.0, 10.0, pk);
     for(k = 0; k < 10; k++)
	cout << breakdown[k] << " "; 
     cout  << endl;
  }*/
return 0;	
}

#include<math.h>
#include<iostream>
#include "../distributionlist.h"
#include "../helperfunctions.h"
#include<iomanip>

using namespace std;
int main()
{
// DistributionList dists;
 double i = 1.0;
 double prev_i;
 double k,l;
 int j;
 
 for(j = 0; j < 10000; j++)
  {
    if( i < MAX_LOG_P)
      cout << j << " "<< log(i) << endl;
    prev_i = i;
    i = i*10;   
 }


// for(j = 1; j < 100; j++)
//  {
 //   i = i + 9*pow(10,-1*j);  
//    k = dists.LogDensity(2, i, BetaDistAlpha(0.999, 10.0), BetaDistBeta(0.999,10.0));
    //l = dists.LogDensity(2, i, BetaDistAlpha(i, 10.0), BetaDistBeta(i, 10.0));
//    cout << j << " " << setprecision(25) <<  i << " "<< 1.0-i << " " << k << " " << l <<  endl;
// }
return 0;
}

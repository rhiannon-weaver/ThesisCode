#include <vector>
#include <iostream>
using namespace std;
int main()
{
  vector<int*> ints;
  int* x;
  int i;
  x = new int[10];
  
  for(i=0; i < 10; i++)
  {
     x[i] = i;
     ints.push_back(&(x[i]));
  } 

  for(i = 0; i < 10; i++)
     cout << *(ints[i]) << " ";
  cout << endl;

  for(i =0 ; i < 10 ; i++)
     cout << x[i] << " ";
  cout << endl;

  ints.clear();

  for(i = 0; i < 10; i++)
    cout << x[i] << " ";
  cout << endl;

  cout << ints.size() << endl;

  return(0);
}

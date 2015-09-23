#include <iostream>
#include <time.h>
#include <stdio.h>
#include <string>
#include <sstream>
using namespace std;
int main()
{
  time_t T;
  struct tm* timeinfo;
  char fmt[11];
  int a;
  int b = 32;
  stringstream str; 
  str << b;
  string b_str(str.str());
  cout << "b = " << b << " = " << b_str <<endl;

  time(&T);
  timeinfo = localtime(&T);
  strftime(fmt, 11, "%d-%m-%Y", timeinfo);
  string s(fmt, 11);
  //char * c = (string("test -d testNetworkGraph.cpp")).c_str;
  cout << fmt << endl << s << endl;
  
  

  a = system("test -d testdir");
  cout << (a==0) << endl;
  a = system("test -f testNetworkGraph.cpp");
  cout << (a==0) << endl;
  a = system("test -d testNetworkGraph.cpp");
  cout << (a==0) << endl;
 // a = system(c);
  cout << (a == 0) << endl;   
 return(0);

}

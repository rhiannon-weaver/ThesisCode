#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <math.h>
#include <stdarg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "distributionlist.h"

using namespace std;


inline string to_string(int i)
{
  stringstream ss;
  ss << i;
  return ss.str();
}


class dat_obj{
public:
  int * fl;
  int * ip;
  int T; 
  short int beginval;
  short int endval;
  short int timezone;
  short int tstar;
  int row_id;

  dat_obj(int N, int rid, ifstream &spikeinfo, ifstream &flinfo, ifstream &ipinfo)
  {
    /*RID is 0-indexed!*/
    /*For processing many lines in a row*/
    GetData(N, rid, spikeinfo, flinfo, ipinfo);
  }
  
  dat_obj(int N, int rid)
  {
    /*RID is 0-indexed!*/
    /*For obtaining data from just one file stream (uses head/tail command line) */
    int flipline;
    int infoline;

    flipline = rid + 1;
    infoline = rid + 2; /*accounts for titles in file*/

    string startcommand = "head -";
    string endipcommand = " ../matrices/matrix.allipcount.0305-0424.txt | tail -1 > _tmp_ip.txt";
    string endflcommand = " ../matrices/matrix.allflcount.0305-0424.txt | tail -1 > _tmp_fl.txt";
    string endinfocommand = " ../matrices/dat_info.txt | tail -1 > _tmp_info.txt";

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
    
    unlink("_tmp_ip.txt");
    unlink("_tmp_fl.txt");
    unlink("_tmp_info.txt");
    
  }

  void Print(int N=-1)
  {
    int i;
    if((N == -1) || (N > T))
      {
	N = T;
      }
    cout << "*****" << endl;
    cout << "Information for data object corresponding to row " << row_id << ": (" << T << " observations)" <<  endl;
    cout << "Begin End TZ Tstar" << endl << beginval << " " << endval << " " << timezone << " " << tstar << endl;
    cout << "First " << N << " IP counts: "; 
    for(i = 0; i < (N-1); i++)
      cout << " " << ip[i];
    cout << endl;
    cout << endl;
    cout << "First " << N << " Flow counts: ";
    for(i = 0; i < N; i++)
      cout << " " << fl[i];
    cout << endl<< "*****" << endl;
    return;
  }
 
  ~dat_obj()
  {
    delete fl;
    delete ip;
  }
private:
  void GetData(int N, int rid, ifstream &spikeinfo, ifstream &flinfo, ifstream &ipinfo)
  {
    int i;
    fl = new int[N];
    ip = new int[N];
    T = N;   
    row_id = rid;
    spikeinfo >> timezone >> tstar >> beginval >> endval;
    //cout << timezone << " " << tstar << " " << beginval << " " << endval << endl;
    for(i = 0; i < N; i++)
      {
	flinfo >> fl[i];
	ipinfo >> ip[i];

	cout << i << " " << fl[i] << " " << ip[i] << ":";

      }
   
    return;
  }
};


void teststream(ifstream *instr, int N)
{
  int i;
  int j;
  for(i = 0 ; i < N; i++)
    {
      *instr >> j;
      cout << j << endl;
    }
  return;

}

int main()
{

  /*  dat_obj testdat(1224, 13);
  cout << "Testdata created" << endl;
  cout.flush();
  testdat.Print(100);
  */

  DistributionList distributions;
  
  int i;
  string i_string;

  i_string = "12345";
  i = atoi(i_string.c_str());
  cout << i << endl;

  /*
  double values[100];
  double logprob[100];
  distributions.Simulate(0, 100, values, 0.0, 2.0);
  
  for(i = 0; i < 100; i++)
    cout << values[i] << " ";
  cout << endl;
  distributions.Simulate(6, 100, values, 0.0, 2.0, -1.0, 1.0);
  for(i = 0; i < 100; i++)
    cout << values[i] << " ";
  cout << endl;
  
  values[0] = -2.0;
  for(i = 0; i < 100; i++)
    {
      logprob[i] = distributions.LogDensity(6, values[i], 0.0, 2.0, -1.0, 1.0);
      cout << logprob[i] <<  " ";
    }
  cout << endl;
  cout << (logprob[0] == (-1.0*INFINITY)) << endl;
  */
  /*
  ifstream instr("../matrices/matrix.allflcount.0305-0424.txt");
  teststream(&instr, 100);

  short int i;
  double a;
  long double b; 
  double powval;
  double x;
  int j;
  a = 0.9;
  b = 100000000000000000.0;
  i = INFINITY;
  powval = pow(a, i);
  cout << powval << endl;
  cout << powval * b << endl;
  cout << i << endl;
  i = i + 1;
  cout << i << endl;
  powval = pow(a, i);
  cout << powval << endl;
  x = 0;
  cout << log(x) << endl;
  
  string stringname("Test string name ");
  string intstring;
  stringstream converter;
  converter << i;
  intstring = converter.str();
  cout << stringname+intstring<< endl;
  */
  return 0;
}

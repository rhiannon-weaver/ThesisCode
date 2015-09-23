#include <map>
#include <vector>
#include <iostream>

using namespace std;
typedef map<short int, int> imap_t;
typedef map<unsigned int, map<short int, int> > omap_t;

int main(){ 
int i, j;

omap_t omap;
omap_t::iterator outer_loop;

imap_t::iterator inner_loop;

omap[12345][2] = 1;
omap[12345][0] = 10;
omap[12345][1] = 5;
omap[67890][2] = 15;
omap[67890][0] = 100;
omap[67890][1] = 50;
omap[34567] = imap_t();

vector<int> tmp(1,12);
cout << tmp.size() << " " << tmp[0] << endl;


for(outer_loop = omap.begin(); outer_loop != omap.end(); outer_loop++)
	{
	  cout << outer_loop->first << ":";	
	  for(inner_loop = outer_loop->second.begin(); inner_loop != outer_loop->second.end(); inner_loop++)
		{
			cout << "(" << inner_loop->first << ","<< inner_loop->second <<") ";
		}
	  cout << endl;
	}	

outer_loop = omap.find(12345);
cout << "Find 12345: " << (outer_loop == omap.end()) << endl;

outer_loop = omap.find(12346);
cout << "Find 12346: " << (outer_loop == omap.end()) << endl;

return(0);
}

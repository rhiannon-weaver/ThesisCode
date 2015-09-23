#ifndef MAINTENANCEFUNCTIONS
#define MAINTENANCEFUNCTIONS
#include <string>
#include <sstream>
using namespace std; 

/* convert a number to a string */
inline string to_string(int i)
{
  stringstream ss;
  ss << i;
  return ss.str();
}

/* comparison funcitons (various) */
int compare_shortint(const void* a, const void* b);
#endif

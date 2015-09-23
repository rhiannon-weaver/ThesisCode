#include "../maintenance_functions.h"


int compare_shortint(const void* a, const void* b)
{
  int val;
  if( *((short int*)(a)) > *((short int*)(b)))
    val = 1;
  else 
    {
      if( *((short int*)(a)) < *((short int*)(b)))
	val = -1;
      else
	val = 0;
    }
  return(val);
}

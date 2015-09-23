#include "test_id.h"

machine_id_t emptyclass::id = 0;
emptyclass::emptyclass()
{
  uid = getID();
}

void emptyclass::printID()
{
  cout << uid << endl;
}

machine_id_t emptyclass::getID()
{
 return id++;
}


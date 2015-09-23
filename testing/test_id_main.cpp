#include <vector>
#include "test_id.h"

int main()
{

int i;
vector<emptyclass*> instances;

for(i = 0; i < 100; i++)
{
	instances.push_back(new emptyclass());
	instances[i]->printID();
}

return 0;
}

#include <iostream>
using namespace std;

void testref(int a, int b, int& sum)
{
 sum = a+b;
}

int main()
{
  int a = 5;
  int b = 4;
  int s = 2;

  cout << a << " " << b << " " << s << endl;

  testref(a, b, s);

  cout << a << " " << b << " " << s << endl;
  
return 0;
}

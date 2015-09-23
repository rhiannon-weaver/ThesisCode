#include <iostream>
using namespace std;

class tst{
  public:
	tst(int a);
 	tst(tst* cpy);
        int getval();
  private:
    int val;

};

tst::tst(int a)
{ val = a;}

tst::tst(tst* cpy)
{val = cpy->val;}

int tst::getval()
{ return(val);}

int main(){

  tst a(12);
  tst b(&a);
  cout << a.getval() << " " << b.getval() << endl;
 return 0;
}
  



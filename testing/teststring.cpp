#include <string>
#include <iostream>

using namespace std;

void PrintString(const char* strval)
{
  string s;
  s = string(strval) + "Info.txt";
  cout << s.c_str() << endl;
}

int main()
{

 string sname;
 string sname2;
 
 sname = "mergesplit";
 sname2 = sname;
 
 cout << sname << " " << sname2 << endl;
 
 sname= "mergesplit2";

 cout << sname << " " << sname2 << endl;

 //PrintString("WaledacInformedSingletons");
 
 return 0;
}

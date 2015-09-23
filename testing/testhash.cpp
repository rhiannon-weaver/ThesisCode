#include <string>
#include <iostream>
#include <map>
using namespace std;

int main(){

map<string, int > hmap;

hmap["one"] = 1;
hmap["two"] = 2;

cout << hmap[hmap.find(2)] << endl;

cout << hmap["two"] << endl;

}

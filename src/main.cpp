#include <iostream>
#include <stdexcept>
#include <string>
#include "ChemicalEquation.h"
int main()
{
	using namespace std;
	string s;
	while(cout << "equation: ",getline(cin,s))
	{
		try{
			cout << ChemicalEquaion::balance(s) << endl;
		}catch(const char *s){
			cout << s << endl;
		}catch(std::exception &e){
			cout << e.what() << endl;
		}
		
	}
	return 0;
}
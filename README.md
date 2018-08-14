# Description
&emsp;A procedure for balancing chemical equations.
# How to use it
**syntax**: See the `main.cpp` file.
```cpp
#include <iostream>
#include <string>
#include <stdexcept>
#include "ChemicalEquaion.h"
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
```
**input format**
* All the chemical formula, +, =, should be separated by Spaces. Example: `CO + O2 = CO2`. 
* Compound states [like (s) (aq) or (g)] are not required, if you input the states, the program will throw a exception, but you can modify this action by the code.
* To enter an ion specify charge after the compound in curly brackets: {+3} or {3+} or {3}. Example: `Fe{3+} + I{-} = Fe{2+} + I2`.
* You can use `(`, `)`, `[`, `]` in the chemical formula.
# What equation can it support
* Type 1: chemical equation
	Eg  : NaBiO3 + MnSO4 + H2SO4 = Na2SO4 + Bi2(SO4)3 + NaMnO4 + H2O
	ans : 10NaBiO3 + 4MnSO4 + 14H2SO4 = 3Na2SO4 + 5Bi2(SO4)3 + 4NaMnO4 + 14H2O
* Type 2: ionic equation
	Eg  : Fe(CN)6{4-} + MnO4{-} + H{+} = Mn{2+} + Fe{3+} + CO2 + NO3{-} + H2O
	ans : 5Fe(CN)6{4-}+61MnO4{-}+188H{+}=61Mn{2+}+5Fe{3+}+30CO2+30NO3{-}+94H2O
* Type 3: chemical equations with multiple balancing results
	Eg  : Fe{2+} + Cl2 + Br{-} = Fe{3+} + Cl{-} + Br2
	ans : (2x1 -2x2)Fe{2+} + (x1)Cl2 + (2x2)Br{-} = (2x1-2x2)Fe{3+} + (2x1)Cl{-} + (x2)Br2

**unsupported types** : Na2Sx + (3x+1)NaClO + (2x-2)NaOH = xNa2SO4 + (3x+1)NaCl + (x-1)H2O
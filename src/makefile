ChemicalEquationBalancer.exe:Fraction.o ChemicalFormula.o ChemicalEquation.o main.o
	g++ -o ChemicalEquationBalancer.exe Fraction.o ChemicalFormula.o ChemicalEquation.o main.o
Fraction.o:Fraction.cpp Fraction.h
	g++ -c Fraction.cpp
ChemicalFormula.o:ChemicalFormula.cpp ChemicalFormula.h pe.h
	g++ -c ChemicalFormula.cpp
ChemicalEquation.o:ChemicalEquation.cpp ChemicalFormula.h ChemicalEquation.h Fraction.h Matrix.h pe.h
	g++ -c ChemicalEquation.cpp
main.o:main.cpp ChemicalEquation.h
	g++ -c main.cpp
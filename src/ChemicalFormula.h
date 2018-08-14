#include <map>
#include <string>
#include <iostream>
#ifndef CHEMICALFORMULA_H_
#define CHEMICALFORMULA_H_
class ChemicalFormula
{
public:
	ChemicalFormula(){};
	std::map<std::string,int> formula;//element to count
    int charge = 0;
	static ChemicalFormula parser(const std::string &chemicalFormula);
	friend std::ostream& operator<<(std::iostream &os,const ChemicalFormula &foo);
};
#endif
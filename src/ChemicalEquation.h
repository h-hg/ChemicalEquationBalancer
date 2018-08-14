#include <iostream>
#include <string>
#include <vector>
#include <set>
#include "Matrix.h"
#include "Fraction.h"
#include "ChemicalFormula.h"
#ifndef CHEMICALeQUATION_H_
#define CHEMICALeQUATION_H_
class ChemicalEquaion
{
public:
	static std::string balance(const std::string &chemicalEquation);
private:
	static void getStrChs(const std::string &chemicalEquation,
			   			  std::vector<std::string> &reactants,
			   			  std::vector<std::string> &resultants);
	static void getChs(const std::vector<std::string> &lhs,
					   const std::vector<std::string> &rhs,
					   std::vector<ChemicalFormula> &reactants,
					   std::vector<ChemicalFormula> &resultants);
	static std::set<std::string> getElements(const std::vector<ChemicalFormula>&chs);
	static Matrix<Fraction> buildMatrix(const std::vector<ChemicalFormula>&reactants,const std::vector<ChemicalFormula>&resultants);
	template <typename T>
	static std::string toString(const T&coeff,const int &num);
	template <typename T>
	static std::string toString(const T&coeff);
	static std::string getEquation(std::vector<std::string> strReactants,
					    std::vector<std::string> strResultants,
						std::vector<std::string> strCoeffs);
	static void getCoeffs(const Matrix<Fraction> &matrix,std::vector<std::string> &strCoeffs);
};
#endif
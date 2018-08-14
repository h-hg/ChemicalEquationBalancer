#include <sstream>
#include "ChemicalEquation.h"
#include "pe.h"
void ChemicalEquaion::getStrChs(const std::string &chemicalEquation,
			   std::vector<std::string> &reactants,
			   std::vector<std::string> &resultants)
{
	std::istringstream iscin(chemicalEquation);
	std::vector<std::string> strChs[2];//string Chemical Formula s
	std::string strCh;
	bool right = false;
	while(iscin >> strCh)
	{
		if(strCh == "=" || strCh == "->")
			right = true;
		else if(strCh != "+"){
			if(right)
				resultants.push_back(strCh);
			else
				reactants.push_back(strCh);
		}
	}
}
void ChemicalEquaion::getChs(const std::vector<std::string> &lhs,
			const std::vector<std::string> &rhs,
			std::vector<ChemicalFormula> &reactants,
			std::vector<ChemicalFormula> &resultants)
{
	//change string chemical to formula
	for(const std::string & foo :lhs)
		reactants.push_back(ChemicalFormula::parser(foo));
	for(const std::string & foo :rhs)
		resultants.push_back(ChemicalFormula::parser(foo));
}
std::set<std::string> ChemicalEquaion::getElements(const std::vector<ChemicalFormula>&chs)
{
	std::set<std::string> ans;
	for(size_t sub = 0;sub < chs.size();++sub)
		for(auto it = chs[sub].formula.cbegin();it != chs[sub].formula.cend();++it)
			ans.insert(it -> first);
	return ans;
}
Matrix<Fraction> ChemicalEquaion::buildMatrix(const std::vector<ChemicalFormula>&reactants,const std::vector<ChemicalFormula>&resultants)
{
	std::set<std::string> elements = getElements(reactants);
	if(elements != getElements(resultants))
		throw pe("Elements on lhs and elements on rhs are different");
	Matrix<Fraction> matrix(elements.size() + 1,reactants.size() + resultants.size());// + 1 for charge
	//coefficient of resultants
	size_t row = 0 ,col = 0;
    for (auto it_elem = elements.begin(); it_elem != elements.end(); ++it_elem, ++row)
    {
        for (col = 0; col < reactants.size(); ++col)
        {
            auto it = reactants[col].formula.find(*it_elem);
			matrix(row,col) = Fraction(it != reactants[col].formula.end() ? it->second : 0);
        }
        for (col = 0; col < resultants.size(); ++col)
        {
            auto it = resultants[col].formula.find(*it_elem);
			matrix(row, reactants.size() + col) = Fraction(it != resultants[col].formula.end() ? -it->second : 0);
        }
    }
	//coefficient of change
	for(col = 0;col < reactants.size();++col)
		matrix(row,col) = Fraction(reactants[col].charge);
	for(col = 0;col < resultants.size();++col)
		matrix(row,reactants.size() + col) = Fraction(-resultants[col].charge);
	return matrix;
}
template <typename T>
std::string ChemicalEquaion::toString(const T&coeff,const int &num)
{
	if(coeff == 0)
		return std::string();
	if(coeff == 1)
		return "x" + std::to_string(num);
	else if(coeff == -1)
		return "-x" + std::to_string(num);
	return std::to_string(coeff) + "x" + std::to_string(num);
}
template <typename T>
std::string ChemicalEquaion::toString(const T&coeff)
{
	if(coeff == 0 || coeff == 1)
		return std::string();
	else if(coeff == -1)
		return "-";
	return std::to_string(coeff);
}
std::string ChemicalEquaion::getEquation(std::vector<std::string> strReactants,
					    std::vector<std::string> strResultants,
						std::vector<std::string> strCoeffs)
{
	std::string ans = strCoeffs[0] + strReactants[0];
	for(size_t i = 1;i < strReactants.size();++i)
		ans += " + " + strCoeffs[i] + strReactants[i];
	ans += " = " + strCoeffs[strReactants.size()] + strResultants[0];
	for(size_t i = 1;i < strResultants.size();++i)
		ans += " + " + strCoeffs[strReactants.size() + i] + strResultants[i];	
	return ans;
}
void ChemicalEquaion::getCoeffs(const Matrix<Fraction> &matrix,std::vector<std::string> &strCoeffs)
{
	std::vector<long long> varState(matrix.getCol(),0);//the default value is 0. positive value stands for pivot element and it's value is the row(start from 1) of pivot elementnegative value stands for free variable and it's value is the ordinal numerber(start from 1) of the free variable
	for(size_t find = 1,i = 0;i < matrix.getRow() && find;++i)
	{
		find = 0;
		for(size_t j = i;j < matrix.getCol();++j)
			if(matrix(i,j)) {
				find = 1;
				varState[j] = i + 1;
				break;
			}
	}
	//find the column of free variable and the number of free variable
	std::vector<std::vector<int>> freeVars;
	std::vector<int> lcms(matrix.getCol(),1);//record the lcm of denumerator of column vectors of free variables，maybe waste some momery
	for(size_t j = 0;j < matrix.getCol();++j)
		if(!varState[j])
		{
			freeVars.push_back(std::vector<int>(matrix.getRow()));
			auto sub = freeVars.size() - 1;//the subcript of current free variable
			//record the column vectors of free variables
			for(size_t i = 0;i < matrix.getRow();++i)
				lcms[sub] = Fraction::lcm(lcms[sub],matrix(i,j).denominator());
			//gain the lcm of denumerator of column vectors of free variables
			for(size_t i = 0;i < matrix.getRow();++i)
				freeVars[sub][i] =  (-lcms[sub]) / matrix(i,j).denominator() * matrix(i,j).numerator();
			varState[j] = - freeVars.size();
		}
	//analyse the free variable
	strCoeffs.resize(matrix.getCol());
	if(freeVars.size() == 0)
		throw pe("There are zero coefficents");
	else if(freeVars.size() == 1)
	{
		for(size_t i = 0;i < matrix.getCol() - freeVars.size();++i)
			if(freeVars[0][i] <= 0)
				throw pe("There are zero or negative coefficents.");
		for(size_t j = 0;j < matrix.getCol();++j)
			if(varState[j] > 0)
				strCoeffs[j] = toString(freeVars[0][varState[j]-1]);
			else
				strCoeffs[j] = toString(lcms[-varState[j]-1]);
	}
	else{
		for(std::string &str:strCoeffs)
			str += '(';
		for(size_t j = 0;j < matrix.getCol();++j)
			if(varState[j] > 0)
				for(size_t k = 0;k < freeVars.size();++k) {
					if(freeVars[k][varState[j]-1] > 0 && strCoeffs[j].empty())
						strCoeffs[j] += "+";
					strCoeffs[j] += toString(freeVars[k][varState[j]-1],k+1);
				}
			else
				strCoeffs[j] += toString(lcms[-varState[j]-1],-varState[j]);
		for(std::string &str:strCoeffs)
			if(str == "(")
				throw pe("There are zero or negative coefficents2.");
			else
				str += ')';
	}
}
std::string ChemicalEquaion::balance(const std::string &chemicalEquation)//The annotated code are using for test
{
	std::vector<std::string> strReactants,strResultants;
	getStrChs(chemicalEquation,strReactants,strResultants);

/*	std::cout << "strReactants(" << strReactants.size() << "): ";
	for(size_t i = 0;i < strReactants.size();++i)
		std::cout << strReactants[i] << ", ";
	std::cout << "\n";
	std::cout << "strResultants(" << strResultants.size() << "): ";
	for(size_t i = 0;i < strResultants.size();++i)
		std::cout << strResultants[i] << ", ";
	std::cout << '\n';*/

	std::vector<ChemicalFormula> reactants,resultants;
	getChs(strReactants,strResultants,reactants,resultants);
	Matrix<Fraction> matrix = buildMatrix(reactants,resultants);

/*	std::cout << "\t" << strReactants[0] << '\t';
	for(size_t i = 1;i < strReactants.size();++i)
		std::cout << "+\t" << strReactants[i] << "\t ";
	std::cout << "=\t" << strResultants[0] << '\t'; 
	for(size_t i = 1;i < strResultants.size();++i)
		std::cout << "+\t" << strResultants[i] << "\t";
	std::cout << std::endl;
	std::set<std::string> elements = getElements(reactants);
	size_t row = 0;
    for (auto it_elem = elements.begin(); it_elem != elements.end(); ++it_elem, ++row)
    {
		std::cout << (*it_elem) << '\t';
		for(size_t col = 0;col < matrix.getCol();++col)
			std::cout << matrix(row,col) << "\t\t";
		std::cout << '\n';
    }
	std::cout << "e\t";
	for(size_t col = 0;col < matrix.getCol();++col)
			std::cout << matrix(row,col) << "\t\t";
	std::cout << '\n';*/

	matrix.rref(matrix.getCol());

	//std::cout << "rref matrix:\n" << matrix << std::endl;

	std::vector<std::string> strCoeffs;
	getCoeffs(matrix,strCoeffs);
	return getEquation(strReactants,strResultants,strCoeffs);
}
#include "ChemicalFormula.h"
#include "pe.h"
#include <stack>
std::ostream& operator<<(std::ostream &os,const ChemicalFormula&foo)
{
	os << "charge: " << foo.charge << '\t';
	for(auto it = foo.formula.begin();it!=foo.formula.end();++it)
		os << (it->first) << ':' << (it->second) << '\t';
	os << std::endl;
	return os;
}
ChemicalFormula ChemicalFormula::parser(const std::string &chemicalFormula)
{
    ChemicalFormula ch;
	std::stack<char> brackets;
	std::stack< std::map<std::string,int> >  elements;
	for(auto pos = chemicalFormula.cbegin();pos!=chemicalFormula.cend();)
	{
		if(isupper(*pos)){
			std::string elemSybol;
			elemSybol += *pos++;
			if(islower(*pos))
				elemSybol += *pos++;
			int elemCount = 0;
			while(isdigit(*pos))
				elemCount = elemCount*10 + (*pos++ -'0');
			if(!elemCount)	elemCount++;
			( elements.empty()?ch.formula:elements.top() )[elemSybol]  += elemCount;
		}
		else if(*pos == '(' || *pos == '[') {
			brackets.push(*pos);
			elements.push(std::map<std::string,int>());
			++pos;
		}
		else if(*pos == ')' || *pos == ']') {
			if( (*pos == ')' && brackets.top()!='(') || (*pos == ']' && brackets.top()!= '['))
				throw pe("mishmatch bracket");
			brackets.pop();
			int multiple = 0;
			for(++pos;isdigit(*pos);)
				multiple = multiple*10 + (*pos++ - '0');
			if(!multiple)	++multiple;
			auto topFormula = elements.top();
			elements.pop();
			std::map<std::string,int>& insertedFormula = elements.empty()?ch.formula:elements.top();//作用于多层括号的时候，如果q为空，说明最外层，直接加在最外层formula
			for(auto it = topFormula.cbegin();it != topFormula.cend();++it)
				insertedFormula[it->first] += it->second *multiple;
		}
		else if(*pos == '{') {
			++pos;
			while(isdigit(*pos))
				ch.charge = ch.charge * 10 + (*pos++ - '0');
			ch.charge = (*pos == '-') ? (ch.charge == 0 ? -1:-ch.charge) : ( (*pos == '+') ? (ch.charge == 0 ? 1:ch.charge) : throw pe("error form of charge") );
			if(*++pos != '}')
				throw pe("error form of charge");
			++pos;
		}
		else
			throw pe("unknown character of chemical formula.");
	}
	return ch;
}
/* using namespace std;
int main()
{
	string s;
	while(cout << "input: ",getline(cin,s))
	{
		try{
			cout << ChemicalFormula::parser(s) << endl;
		}catch(const char *s){
			cout << s << endl;
		}
	}
	return 0;
} */
#include "Fraction.h"
#include <stdexcept>
std::istream& operator>>(std::istream &is,Fraction &fraction)
{
	char c;
	is >> fraction.nu;
	c = is.peek();
	if(c == '/')
		is.ignore() >> fraction.de;
	else
		fraction.de = 1;
	if(fraction.de == 0)
		throw std::domain_error("denominator can't be 0");
    fraction.simplify();
	return is;
}
std::ostream& operator<<(std::ostream &os,const Fraction &fraction)
{
 	os << fraction.nu;
	if(fraction.nu != 0 && fraction.de != 1)
		os << '/' << fraction.de;
	return os;
}
void Fraction::assign(long long& newNu, long long& newDe)
{
    long long g = gcd(newNu > 0 ? newNu : -newNu, newDe);
    newNu /= g;
    newDe /= g;
    nu = static_cast<int>(newNu);
    de = static_cast<int>(newDe);
}
void Fraction::simplify()
{
    if(!de) throw std::domain_error("denominator can't be 0");
    if (de < 0) { nu = -nu; de = -de;}
    int g = gcd(nu > 0 ? nu : -nu, de);
    nu /= g;
    de /= g;
}
Fraction::Fraction(int numerator,int denominator):nu(numerator),de(denominator)
{
    simplify();
}
Fraction::Fraction(const Fraction&fraction):nu(fraction.nu),de(fraction.de){}
Fraction& Fraction::operator+=(const Fraction&rhs)
{
    long long newDe = static_cast<long long>(de) * rhs.de;
    long long newNu = static_cast<long long>(nu) * rhs.de + static_cast<long long>(rhs.nu) * de;
    assign(newNu, newDe);
    if (nu != newNu || de != newDe) throw std::overflow_error("operator+ overflows");
    return *this;
}
Fraction& Fraction::operator-=(const Fraction&rhs)
{
    long long newDe = static_cast<long long>(de) * rhs.de;
    long long newNu = static_cast<long long>(nu) * rhs.de - static_cast<long long>(rhs.nu) * de;
    assign(newNu, newDe);
    if (nu != newNu || de != newDe) throw std::overflow_error("operator- overflows");
    return *this;
}
Fraction& Fraction::operator*=(const Fraction&rhs)
{
    long long newDe = static_cast<long long>(de) * rhs.de;
    long long newNu = static_cast<long long>(nu) * rhs.nu;
    assign(newNu, newDe);
    if (nu != newNu || de != newDe) throw std::overflow_error("operator* overflows");
    return *this;
}
Fraction& Fraction::operator/=(const Fraction&rhs)
{
    if (rhs.nu == 0) throw std::domain_error("Divide by 0");
    long long newDe = static_cast<long long>(de) * rhs.nu;
    long long newNu = static_cast<long long>(nu) * rhs.de;
    if (newDe < 0) { newDe = -newDe; newNu = -newNu; }
    assign(newNu, newDe);
    if (nu != newNu || de != newDe) throw std::overflow_error("operator* overflows");
    return *this;
}
Fraction& Fraction::operator=(const Fraction&rhs)
{
    nu = rhs.nu;
    de = rhs.de;
    return *this;
}
Fraction Fraction::operator+(const Fraction&rhs)const
{
    Fraction lhs = *this;
    return lhs += rhs;
}
Fraction Fraction::operator-(const Fraction&rhs)const
{
    Fraction lhs = *this;
    return lhs -= rhs;
}
Fraction Fraction::operator*(const Fraction&rhs)const
{
    Fraction lhs = *this;
    return lhs *= rhs;
}
//Fraction operator%(const Fraction&rhs)const;
Fraction Fraction::operator/(const Fraction&rhs)const
{
    Fraction lhs = *this;
    return lhs /= rhs;
}
Fraction Fraction::operator-()const
{
    return Fraction(-nu,de);
}
bool Fraction::operator==(const Fraction&rhs)const
{
    return (nu == rhs.nu && de == rhs.de);
}
bool Fraction::operator!=(const Fraction&rhs)const
{
    return !(nu == rhs.nu && de == rhs.de);
}
bool Fraction::operator<(const Fraction&rhs)const
{
    return static_cast<long long>(nu) * rhs.de < static_cast<long long>(de) * rhs.nu;
}
bool Fraction::operator>(const Fraction&rhs)const
{
    return static_cast<long long>(nu) * rhs.de > static_cast<long long>(de) * rhs.nu;
}
bool Fraction::operator<=(const Fraction&rhs)const
{
    return static_cast<long long>(nu) * rhs.de <= static_cast<long long>(de) * rhs.nu;
}
bool Fraction::operator>=(const Fraction&rhs)const
{
    return static_cast<long long>(nu) * rhs.de >= static_cast<long long>(de) * rhs.nu;
}
bool Fraction::operator!()const
{
    return !nu;
}
Fraction::operator bool()const
{
    return nu;
}
Fraction::operator int()const
{
    return nu/de;
}
Fraction::operator float()const
{
    return static_cast<float>(nu)/de;
}
Fraction::operator double()const
{
    return static_cast<double>(nu)/de;
}
Fraction::operator long double()const
{
    return static_cast<long double>(nu)/de;
}
int Fraction::numerator()const
{
    return nu;
}
int Fraction::denominator()const
{
    return de;
}
void Fraction::reset(int numerator,int denominator)
{
    nu = numerator;
    de = denominator;
    simplify();
}
Fraction Fraction::abs()const
{
    return Fraction(nu > 0 ? nu :-nu,de);
}
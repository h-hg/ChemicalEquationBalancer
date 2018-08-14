#include <iostream>
#ifndef FRACTION_H_
#define FRACTION_H_
class Fraction
{
    friend std::ostream& operator<<(std::ostream &os,const Fraction &fraction);
	friend std::istream& operator>>(std::istream &is,Fraction &fraction);
public:
	Fraction(int numerator = 0,int denominator = 1);
    Fraction(const Fraction&fraction);

    Fraction& operator=(const Fraction&rhs);
	Fraction& operator+=(const Fraction&rhs);
	Fraction& operator-=(const Fraction&rhs);
	Fraction& operator*=(const Fraction&rhs);
	Fraction& operator/=(const Fraction&rhs);
	Fraction operator+(const Fraction&rhs)const;
	Fraction operator-(const Fraction&rhs)const;
	//Fraction operator%(const Fraction&rhs)const;
	Fraction operator/(const Fraction&rhs)const;
    Fraction operator*(const Fraction&rhs)const;
    Fraction operator-()const;

    bool operator==(const Fraction&rhs)const;
    bool operator!=(const Fraction&rhs)const;
    bool operator<(const Fraction&rhs)const;
    bool operator>(const Fraction&rhs)const;
    bool operator<=(const Fraction&rhs)const;
    bool operator>=(const Fraction&rhs)const;
    bool operator!()const;

    operator bool()const;
    operator int()const;
    operator float()const;
    operator double()const;
    operator long double()const;
    
    int numerator()const;
    int denominator()const;
    void reset(int numerator = 0,int denominator = 1);
    Fraction abs()const;
    template<typename T>
    static T gcd(const T&,const T&);
    template<typename T>
    static T lcm(const T&,const T&);
private:
    //static long long gcd(long long&,long long&);
    void simplify();
    void assign(long long& newNu, long long& newDe);
    int nu;
    int de;
};
template<typename T>
T Fraction::gcd(const T &t1,const T &t2)//The binary GCD algorithm, also known as Stein's algorithm.
{
    T x = t1,y = t2;
    if(!x)  return y;
    if(!y)  return x;
    int times = 0;
    while(!(x&1) && !(y&1))
    {
        ++times;
        x >>= 1;
        y >>= 1;
    }
    while(!(x&1))   x>>=1;
    while(!(y&1))   y>>=1;
    while(x != y)
        if(x > y){
            x -= y;
            while(!(x&1))   x>>=1;
        }
        else{
            y -= x;
            while(!(y&1))   y>>=1;
        }
    return x<<times;
}
template<typename T>
T Fraction::lcm(const T &t1,const T &t2)
{
    return t1 / gcd(t1,t2) * t2;
}
#endif
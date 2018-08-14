#include <iostream>
#include <cstring>
#ifndef Matrix_H_
#define Matrix_H_
template <typename T>
class Matrix
{
	template <typename U>
	friend std::ostream& operator<<(std::ostream &os,const Matrix<U> &matrix);
	template <typename U>
	friend std::istream& operator>>(std::istream &is,Matrix<U> &matrix);
private:
	void build();
	void clear();
	void copy(const Matrix<T> &matrix);
	void movecopy(Matrix<T> &matrix);

	size_t row;
	size_t column;
	T **data;
public:
	Matrix(size_t row,size_t column);
	Matrix(const Matrix<T> &matrix);
	Matrix(Matrix<T> &&matrix);
	~Matrix();
	Matrix& operator=(const Matrix<T> &matrix);
	Matrix& operator=(Matrix<T> &&matrix);

	T& operator()(size_t row,size_t column);
	T operator()(size_t row,size_t column)const;
	size_t getCol()const{return column;}
	size_t getRow()const{return row;}

	void rowInterchange(const size_t &row1,const size_t &row2);
	void rref(size_t limitCol);
};
#endif

template <typename U>
std::ostream& operator<<(std::ostream &os,const Matrix<U> &matrix)
{
	for(size_t i = 0;i < matrix.row;++i)
	{
		for(size_t j = 0; j < matrix.column - 1;++j)
			os << matrix.data[i][j] << '\t';
		os << matrix.data[i][matrix.column-1] << '\n';
	}
	return os;
}
template <typename U>
std::istream& operator>>(std::istream &is,Matrix<U> &matrix)
{
	for(size_t i = 0;i < matrix.row;++i)
		for(size_t j = 0;j < matrix.column;++j)
			is >> matrix.data[i][j];
	return is;
}
template<typename T>
void Matrix<T>::build()
{
	data = new T*[row];
	for(size_t sub = 0;sub < row;++sub)
		data[sub] = new T[column];
}
template<typename T>
void Matrix<T>::clear()
{
	for(size_t sub = 0;sub < row;++sub)
		delete[] data[sub];
	delete[] data;
}
template<typename T>
void Matrix<T>::copy(const Matrix<T> &matrix)
{
	row = matrix.row;
	column = matrix.column;
	build();
	for(auto sub = 0;sub < row;++sub)
		memcpy(data[sub],matrix.data[sub],sizeof(T)*column);
}
template<typename T>
void Matrix<T>::movecopy(Matrix<T> &matrix)
{
	row = matrix.row;
	column = matrix.column;
	data = matrix.data;
	matrix.data = nullptr;
}

template<typename T>
Matrix<T>::Matrix(size_t row,size_t column):row(row),column(column)
{
	build();
}
template<typename T>
Matrix<T>::Matrix(const Matrix<T> &matrix)
{
	copy(matrix);
}
template<typename T>
Matrix<T>::Matrix(Matrix<T> &&matrix)
{
	movecopy(matrix);
}
template<typename T>
Matrix<T>::~Matrix()
{
	clear();
}
template <typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T> &matrix)
{
	if(this == &matrix)
		return *this;
	clear();
	copy(matrix);
	return *this;
}
template <typename T>
Matrix<T>& Matrix<T>::operator=(Matrix<T> &&matrix)
{
	if(this == &matrix)
		return *this;
	movecopy(matrix);
	return *this;
}
template <typename T>
T& Matrix<T>::operator()(size_t row,size_t column)
{
	return data[row][column];
}
template <typename T>
T Matrix<T>::operator()(size_t row,size_t column)const
{
	return data[row][column];
}
template<typename T>
void Matrix<T>::rowInterchange(const size_t &row1,const size_t &row2)
{
	auto temp = data[row1];
	data[row1] = data[row2];
	data[row2] = temp;
}

template <typename T>
void Matrix<T>::rref(size_t limitCol)
{
	size_t currentRow = 0,currentCol = 0;
	//the followed steps is converting the matrix to Row-Echelon Form
	for(;currentRow < row && currentCol < limitCol;++currentCol)
	{
		size_t non_zero_row = currentRow;
		//find first row which has non-zero entry at column currentCol
		for(;non_zero_row < row;++non_zero_row)
			if(data[non_zero_row][currentCol])
				break;
		if(non_zero_row == row)//can't find the non-zero iterm, enter the next column
			continue;
		if(non_zero_row != currentRow)
			rowInterchange(non_zero_row,currentRow);
		//minus
		for(auto i = currentRow + 1;i < row;++i)
		{
			auto foo = data[i][currentCol] / data[currentRow][currentCol];
			for(auto j = currentCol;j < column;++j)
				data[i][j] -=  foo * data[currentRow][j];
		}
		++currentRow;
	}
	//std::cout << (*this) << std::endl;
	//the followed steps is converting the matrix to Reduced row echelon form
  	while(currentRow--)
	{
		//find first non-zero entry 
		for(currentCol = 0;currentCol < limitCol;++currentCol)
			if(data[currentRow][currentCol])
				break;
		//minus
		for(size_t i = 0; i < currentRow;++i)
		{
			auto foo = data[i][currentCol] / data[currentRow][currentCol];
			//std::cout << "foo: " << foo << std::endl;
			if(data[i][currentCol])
				for(size_t j = currentCol;j < column;++j)
					data[i][j] -= foo * data[currentRow][j];
		}
		//std::cout << (*this) << std::endl;
		// unit 1
		auto foo = data[currentRow][currentCol];
		for(size_t j = currentCol;j < column;++j)
			data[currentRow][j] /= foo;
		//std::cout << (*this) << std::endl;
	} 
}
/*#include "Fraction.cpp"
using namespace std;
int main()
{
	try{
	freopen("in.txt","r",stdin);
	freopen("out_new.txt","w",stdout);
	size_t row,column;
	cin >> row >> column;
	Matrix<Fraction> *p = new Matrix<Fraction>(row,column);
	cin >> (*p);
	cout << (*p) << endl;
	p->rref(p->column);
	cout << (*p);
	}catch(...){
		cout << "exception";
	}
	return 0;
}*/
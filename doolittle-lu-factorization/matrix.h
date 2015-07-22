/*
 * matrix.h
 *
 *  Created on: Jul 16, 2015
 *      Author: nino
 */

#ifndef MATRIX_H_
#define MATRIX_H_



#include <fstream>
#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <fcntl.h>
#include <string.h>
#include <ctype.h>

#define LOOP_UNROLL 10

using namespace std;

template <class T>
inline void SWAP(T &a, T &b)
{
	T tmp = a;
	a = b;
	b = tmp;
}

template <class T>
inline const T &MAX(const T &a, const T &b)
{
	return (a >= b ? a : b);
}

template <class T>
inline const T &MIN(const T &a, const T &b)
{
	return (a <= b ? a : b);
}

template <class T>
class Matrix {
	private:
		int _n;
		int _m;
		T **v;
	public:
		Matrix();
		explicit Matrix(const int n, const int m);
		Matrix(const int n, const int m, const T &a);
		Matrix(const int n, const int m, const T **a);
		Matrix(const Matrix &rhs);
		Matrix & operator=(const Matrix &rhs);
		inline T* operator[](const int i) const;
		T** getV() const;
		inline const int n_rows() const;
		inline const int n_cols() const;
		inline void swapRows(const int i1, const int i2);
		inline void swapRowsNaive(const int i1, const int i2);
		inline void swapCols(const int j1, const int j2);
		inline void swapColsNaive(const int j1, const int j2);
		inline void multiplyRow(const int i, const T &c);
		inline void multiplyRowNaive(const int i, const T &c);
		inline void multiplyCol(const int j, const T &c);
		inline void multiplyColNaive(const int j, const T &c);
		inline void addRows(const T &c, const int i1, const int i2);
		void print() const;
		~Matrix();
};


template <class T>
Matrix<T>::Matrix() : _n(0), _m(0), v(NULL) {}

template <class T>
Matrix<T>::Matrix(const int n, const int m) : _n(n), _m(m), v(n > 0 ? new T*[n] : NULL)
{
	// the matrix is stored in contiguous memory locations
	int n_elems = n * m;
	if(v)
		v[0] = n_elems > 0 ? new T[n_elems] : NULL;
	for(int i = 1; i < n; ++i)
		v[i] = v[i - 1] + m;
}

template <class T>
Matrix<T>::Matrix(const int n, const int m, const T &a) : _n(n), _m(m), v(n > 0 ? new T*[n] : NULL)
{
	// the matrix is stored in contiguous memory locations
	int n_elems = n * m;
	if(v)
		v[0] = n_elems > 0 ? new T[n_elems] : NULL;
	for(int i = 1; i < n; ++i)
		v[i] = v[i - 1] + m;
	for(int i = 0; i < n; ++i)
		for(int j = 0; j < m; ++j)
			v[i][j] = a;
}

template <class T>
Matrix<T>::Matrix(const int n, const int m, const T **a) : _n(n), _m(m), v(n > 0 ? new T*[n] : NULL)
{
	int n_elems = n * m;
	if(v)
		v[0] = n_elems > 0 ? new T[n_elems] : NULL;
	for(int i = 1; i < n; ++i)
		v[i] = v[i - 1] + m;
	for(int i = 0; i < n; ++i)
		for(int j = 0; j < m; ++j)
			v[i][j] = a[i][j];
}

template <class T>
Matrix<T>::Matrix(const Matrix &rhs) : _n(rhs._n), _m(rhs._m), v(rhs._n > 0 ? new T*[rhs._n] : NULL)
{
	int n_elems = _n * _m;
	if(v)
		v[0] = n_elems > 0 ? new T[n_elems] : NULL;
	for(int i = 1; i < _n; i++)
		v[i] = v[i - 1] + _m;
	for(int i = 0; i < _n; i++)
		for(int j = 0; j < _m; j++)
			v[i][j] = rhs[i][j];
}

template <class T>
Matrix<T>::~Matrix()
{
	if(v != NULL)
	{
		delete[] v[0];
		delete[] v;
	}
}

template <class T>
Matrix<T> & Matrix<T>::operator=(const Matrix<T> &rhs)
{
	if(this != &rhs)
	{
		if(_n != rhs._n || _m != rhs._m)
		{
			if(v != NULL)
			{
				delete [] v[0];
				delete [] v;
			}
			_n = rhs._n;
			_m = rhs._m;
			v = _n > 0 ? new T*[_n] : NULL;
			int n_elems = _n * _m;
			if(v)
				v[0] = n_elems > 0 ? new T[n_elems] : NULL;
			for(int i = 1; i < _n; i++)
				v[i] = v[i - 1] + _m;
		}
		for(int i = 0; i < _n; ++i)
			for(int j = 0; j < _m; j++)
				v[i][j] = rhs[i][j];
	}
	return *this;
}

template <class T>
inline T* Matrix<T>::operator[](const int i) const
{
	// no check for index bounds to minimize branching
	return v[i];
}

template <class T>
T** Matrix<T>::getV() const
{
	return v;
}


template <class T>
inline const int Matrix<T>::n_rows() const
{
	return _n;
}

template <class T>
inline const int Matrix<T>::n_cols() const
{
	return _m;
}

template <class T>
inline void Matrix<T>::swapRows(const int i1, const int i2)
{
	int k = 0;
	for( ; k < _m - LOOP_UNROLL + 1; k += LOOP_UNROLL)
	{
		SWAP(v[i1][k], v[i2][k]);
		SWAP(v[i1][k + 1], v[i2][k + 1]);
		SWAP(v[i1][k + 2], v[i2][k + 2]);
		SWAP(v[i1][k + 3], v[i2][k + 3]);
		SWAP(v[i1][k + 4], v[i2][k + 4]);
		SWAP(v[i1][k + 5], v[i2][k + 5]);
		SWAP(v[i1][k + 6], v[i2][k + 6]);
		SWAP(v[i1][k + 7], v[i2][k + 7]);
		SWAP(v[i1][k + 8], v[i2][k + 8]);
		SWAP(v[i1][k + 9], v[i2][k + 9]);
	}
	for( ; k < _m; ++k)
		SWAP(v[i1][k], v[i2][k]);
}

template <class T>
inline void Matrix<T>::swapRowsNaive(const int i1, const int i2)
{
	for(int k = 0; k < _m; k++)
		SWAP(v[i1][k], v[i2][k]);
}

template <class T>
inline void Matrix<T>::swapCols(const int j1, const int j2)
{
	int k = 0;
	for( ; k < _n - LOOP_UNROLL + 1; k += LOOP_UNROLL)
	{
		SWAP(v[k][j1], v[k][j2]);
		SWAP(v[k + 1][j1], v[k + 1][j2]);
		SWAP(v[k + 2][j1], v[k + 2][j2]);
		SWAP(v[k + 3][j1], v[k + 3][j2]);
		SWAP(v[k + 4][j1], v[k + 4][j2]);
		SWAP(v[k + 5][j1], v[k + 5][j2]);
		SWAP(v[k + 6][j1], v[k + 6][j2]);
		SWAP(v[k + 7][j1], v[k + 7][j2]);
		SWAP(v[k + 8][j1], v[k + 8][j2]);
		SWAP(v[k + 9][j1], v[k + 9][j2]);
	}
	for( ; k < _n; ++k)
		SWAP(v[k][j1], v[k][j2]);
}

template <class T>
inline void Matrix<T>::swapColsNaive(const int j1, const int j2)
{
	for(int k = 0; k < _n; ++k)
			SWAP(v[k][j1], v[k][j2]);

}

template <class T>
inline void Matrix<T>::multiplyRow(const int i, const T &c)
{
	int k = 0;
	for( ; k < _m - LOOP_UNROLL + 1; k += LOOP_UNROLL)
	{
		v[i][k] *= c;
		v[i][k + 1] *= c;
		v[i][k + 2] *= c;
		v[i][k + 3] *= c;
		v[i][k + 4] *= c;
		v[i][k + 5] *= c;
		v[i][k + 6] *= c;
		v[i][k + 7] *= c;
		v[i][k + 8] *= c;
		v[i][k + 9] *= c;
	}
	for( ; k < _m; ++k)
		v[i][k] *= c;
}

template <class T>
inline void Matrix<T>::multiplyRowNaive(const int i, const T &c)
{
	for(int k = 0; k < _m; k++)
		v[i][k] *= c;
}

template <class T>
inline void Matrix<T>::multiplyCol(const int j, const T &c)
{
	int k = 0;
	for( ; k < _n - LOOP_UNROLL + 1; k += LOOP_UNROLL)
	{
		v[k][j] *= c;
		v[k + 1][j] *= c;
		v[k + 2][j] *= c;
		v[k + 3][j] *= c;
		v[k + 4][j] *= c;
		v[k + 5][j] *= c;
		v[k + 6][j] *= c;
		v[k + 7][j] *= c;
		v[k + 8][j] *= c;
		v[k + 9][j] *= c;
	}
	for( ; k < _n; k++)
		v[k][j] *= c;
}

template <class T>
inline void Matrix<T>::multiplyColNaive(const int j, const T &c)
{
	for(int k = 0; k < _n; k++)
		v[k][j] *= c;
}

template <class T>
inline void Matrix<T>::addRows(const T &c, const int i1, const int i2)
{
	int j = 0;
	for(j = 0; j < _m - LOOP_UNROLL + 1; j+= LOOP_UNROLL)
	{
		v[i2][j] += v[i1][j] * c;
		v[i2][j + 1] += v[i1][j + 1] * c;
		v[i2][j + 2] += v[i1][j + 2] * c;
		v[i2][j + 3] += v[i1][j + 3] * c;
		v[i2][j + 4] += v[i1][j + 4] * c;
		v[i2][j + 5] += v[i1][j + 5] * c;
		v[i2][j + 6] += v[i1][j + 6] * c;
		v[i2][j + 7] += v[i1][j + 7] * c;
		v[i2][j + 8] += v[i1][j + 8] * c;
		v[i2][j + 9] += v[i1][j + 9] * c;
	}
	for( ; j < _m; j++)
		v[i2][j] += v[i1][j] * c;
}

template <class T>
void Matrix<T>::print() const
{
	for(int i = 0; i < _n; ++i)
	{
		for(int j = 0; j < _m; ++j)
			printf("\t%.5lf", v[i][j]);
		cout << endl;
	}
}


template <class T>
class Vector {
	private:
		int _n;
		T *v;
	public:
		Vector();
		explicit Vector(const int n);
		Vector(const int n, const T &a);
		Vector(const int n, const T *a);
		Vector(const Vector &rhs);
		Vector & operator=(const Vector &rhs);
		inline T & operator[](const int i) const;
		inline const int size() const;
		inline void swapElements(const int i1, const int i2);
		inline void multiply(const T &c);
		void print() const;
		~Vector();
};


template <class T>
Vector<T>::Vector() : _n(0), v(NULL) {}

template <class T>
Vector<T>::Vector(const int n) : _n(n), v(n > 0 ? new T[n] : NULL) {}

template <class T>
Vector<T>::Vector(const int n, const T &a) : _n(n), v(n > 0 ? new T[n] : NULL)
{
	for(int i = 0; i < n; ++i)
		v[i] = a;
}

template <class T>
Vector<T>::Vector(const int n, const T *a) : _n(n), v(n > 0 ? new T[n] : NULL)
{
	for(int i = 0; i < n; ++i)
		v[i] = *a++;
}

template <class T>
Vector<T>::Vector(const Vector<T> &rhs) : _n(rhs._n), v(rhs._n > 0 ? new T[rhs._n] : NULL)
{
	for(int i = 0; i < rhs._n; ++i)
		v[i] = rhs[i];
}

template <class T>
Vector<T> & Vector<T>::operator=(const Vector<T> &rhs)
{
	if(this != &rhs)
	{
		if(v != NULL)
			delete [] v;
		if(_n != rhs._n)
		{
			_n = rhs._n;
			v = _n > 0 ? new T[rhs._n] : NULL;
		}
		for(int i = 0; i < _n; ++i)
			v[i] = rhs[i];
	}
	return *this;
}

template <class T>
inline T & Vector<T>::operator[](const int i) const
{
	return v[i];
}

template <class T>
inline const int Vector<T>::size() const
{
	return _n;
}

template <class T>
inline void Vector<T>::swapElements(const int i1, const int i2)
{
	SWAP(v[i1], v[i2]);
}

template <class T>
inline void Vector<T>::multiply(const T &c)
{
	int i = 0;
	for( ; i < _n - LOOP_UNROLL + 1; i += LOOP_UNROLL)
	{
		v[i] *= c;
		v[i+1] *= c;
		v[i+2] *= c;
		v[i+3] *= c;
	}
	for( ; i < _n; ++i)
		v[i] *= c;
}

template <class T>
void Vector<T>::print() const
{
	for(int i =0 ; i < _n; ++i)
		cout << v[i] << " ";
	cout << endl;
}

template <class T>
Vector<T>::~Vector()
{
	if(v != NULL)
	delete[] (v);
}

/* Some useful typedefs to keep things clear */
typedef Matrix<double> MatDouble;
typedef Matrix<float> MatFloat;
typedef Matrix<int> MatInt;
typedef Matrix<long> MatLong;
typedef Vector<double> VecDouble;
typedef Vector<float> VecFloat;
typedef Vector<int> VecInt;
typedef Vector<long> VecLong;



#endif /* MATRIX_H_ */

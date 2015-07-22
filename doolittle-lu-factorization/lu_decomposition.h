#ifndef LU_DECOMPOSITION_H_
#define LU_DECOMPOSITION_H_
#include <cmath>

#include "matrix.h"

void doolittle(MatDouble &a, MatDouble &L, MatDouble &U)
{
	int n = a.n_rows();
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
		{
			L[i][j] = 0.0;
			U[i][j] = 0.0;
		}

	for (int k = 0; k < n; k++)
	{
		L[k][k] = 1.0;
		for (int j = k; j < n; j++)
		{
			double s = 0.0;
			for (int m = 0; m < k; m++)
			{
				s += L[k][m] * U[m][j];
			}
			U[k][j] = a[k][j] - s;
		}
		for (int i = k + 1; i < n; i++)
		{
			double s = 0.0;
			for (int m = 0; m < k; m++)
			{
				s += L[i][m] * U[m][k];
			}
			L[i][k] = (a[i][k] - s) / U[k][k];
		}
	}
}

void solveForward(MatDouble &L, VecDouble &b, VecDouble &y)
{
	int n = L.n_rows();
	for (int i = 0; i < n; i++)
		y[i] = 0.0;

	y[n - 1] = b[n - 1] / L[n - 1][n - 1];

	for (int i = 0; i < n - 1; i++)
	{
		double s = 0.0;
		for (int j = 0; j < i; j++)
		{
			s += L[i][j] * y[j];
		}
		y[i] = (b[i] - s) / L[i][i];
	}
}

void solveBackward(MatDouble &U, VecDouble &y, VecDouble &x)
{
	int n = U.n_rows();
	for (int i = 0; i < n; i++)
		x[i] = 0.0;

	x[n - 1] = y[n - 1] / U[n - 1][n - 1];

	for (int i = n - 2; i >= 0; i--)
	{
		double s = 0.0;
		for (int j = i + 1; j < n; j++)
		{
			s += U[i][j] * x[j];
		}
		x[i] = (y[i] - s) / U[i][i];
	}
}

void solve(MatDouble &L, MatDouble &U, VecDouble &b, VecDouble &x)
{
	int n = b.size();
	VecDouble y(n);
	solveForward(L, b, y);
	solveBackward(U, y, x);
}

#endif

/*
 * main.cpp
 *
 *  Created on: Jul 13, 2015
 *      Author: nino
 */

//#include "../gaussj.h"
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
#include<cpucounters.h>
using namespace std;
SystemCounterState before_state, after_state;

void gaussj(const int n, const int m, double **a, double **b)
{
	int ipiv, i, j, k, l;
	double factor;
	double pivotmax;
	double tmp;
	for (i = 0; i < n; i++)
	{
		ipiv = i;
		for (j = i + 1; j < n; j++)
		{
			if (fabs(a[j][i]) > fabs(a[ipiv][i]))
			{
				ipiv = j;
			}
		}

		if (a[ipiv][ipiv] == 0)
			throw("Singular");

		if (ipiv != i)
		{
			#pragma ivdep
			for (j = 0; j < n; j++)
			{
				tmp = a[ipiv][j];
				a[ipiv][j] = a[i][j];
				a[i][j] = tmp;
			}

			#pragma ivdep
			for (l = 0; l < m; l++)
			{
				tmp = b[ipiv][l];
				b[ipiv][l] = b[i][l];
				b[i][l] = tmp;
			}
		}

		for (j = 0; j < i; j++)
		{
			factor = -a[j][i] / a[i][i];

			#pragma ivdep
			for (k = j + 1; k < n; k++)
			{
				a[j][k] += a[i][k] * factor;
			}

			#pragma ivdep
			for (k = 0; k < m; k++)
			{
				b[j][k] += b[i][k] * factor;
			}

		}

		for (j = i + 1; j < n; j++)
		{
			factor = -a[j][i] / a[i][i];

			#pragma ivdep
			for (k = 0; k < n; k++)
			{
				a[j][k] += a[i][k] * factor;
			}

			#pragma ivdep
			for (k = 0; k < m; k++)
			{
				b[j][k] += b[i][k] * factor;
			}
		}
	}

	for (i = 0; i < n; i++)
	{
		#pragma ivdep
		for (j = 0; j < m; j++)
			b[i][j] /= a[i][i];
		a[i][i] = 1.0;
	}
}

int main(int argc, char* argv[])
{
	if (argc <= 1)
	{
		cerr << "Error: No arguments provided!" << endl;
		exit(1);
	}

	const int N = atoi(argv[1]);
	const int r = atoi(argv[2]);

	PCM *m = PCM::getInstance();
	PCM::ErrorCode returnResult = m->program();
	if (returnResult != PCM::Success)
	{
		std::cerr << "Intel's PCM couldn't start" << std::endl;
		std::cerr << "Error code: " << returnResult << std::endl;
		exit(1);
	}

	double **a = new double*[N];
	a[0] = new double[N * N];
	for (int i = 1; i < N; i++)
		a[i] = a[i - 1] + N;

	double **b = new double*[N];
	b[0] = new double[N * r];
	for (int i = 1; i < N; i++)
		b[i] = b[i - 1] + r;

	srand(time(NULL));
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
			a[i][j] = rand() % 100;
		for (int j = 0; j < r; j++)
			b[i][j] = rand() % 40;
	}

	clock_t begin, end;
	SystemCounterState before_state, after_state;
	try
	{
		before_state = getSystemCounterState();
		begin = clock();
		gaussj(N, r, a, b);
		end = clock();
		after_state = getSystemCounterState();

	} catch (const char* e)
	{
		cout << e << endl;
	}
	double time = (double)(end - begin);

	cout << endl << "-------------------------------------------------------------------------" << endl;
	cout << "Total number of retired instructions: " << getInstructionsRetired(before_state, after_state) << endl;
	cout << "Consumed energy by CPU in mW: " << getConsumedJoules(before_state, after_state)/(time*1000.0) << endl;
	cout << "Consumed energy by DRAM in mW: " << getDRAMConsumedJoules(before_state, after_state)/(time*1000.0) << endl;
	cout << "Bytes read from DRAM controller: " << getBytesReadFromMC(before_state, after_state) << endl;
	cout << "Bytes written to DRAM controller: " << getBytesWrittenToMC(before_state, after_state) << endl;
	cout << "Number of CPU cycles: " << getCycles(before_state, after_state) << endl;
	cout << "Instructions per cycle (IPC): " << getIPC(before_state, after_state) << endl;
	cout << "Maximum possible IPC: " << m->getMaxIPC() << endl;
	cout << "Average number of retired instructions per time interval: " << getTotalExecUsage(before_state, after_state) << endl;
	cout << "Average core frequency: " << getAverageFrequency(before_state, after_state) / 1000000000 << " GHz" <<  endl;
	cout << "Cycles lost due to L2 cache misses but still hitting L3 cache: " << getCyclesLostDueL2CacheMisses(before_state, after_state) << endl;
	cout << "Cycles lost due to L3 cache misses: " << getCyclesLostDueL3CacheMisses(before_state, after_state) << endl;
	cout << "L2 cache misses: " << getL2CacheMisses(before_state, after_state) << endl;
	cout << "L3 cache misses: " << getL3CacheMisses(before_state, after_state) << endl;
	cout << "L2 cache hit rate: " << getL2CacheHitRatio(before_state, after_state) << endl;
	cout << "L3 cache hit rate: " << getL3CacheHitRatio(before_state, after_state) << endl;
	cout << "-------------------------------------------------------------------------" << endl;

	cout << N << " " << time * 1000.0 / CLOCKS_PER_SEC << endl;

	m->cleanup();

	delete[] a[0];
	delete[] a;
	delete[] b[0];
	delete[] b;

	return 0;
}

/*
 * main.cpp
 *
 *  Created on: Jul 13, 2015
 *      Author: nino
 */

#include "../lu_decomposition.h"

#include<iostream>
#include<ctime>
#include<cstdlib>
#include<cpucounters.h>
using namespace std;

int main(int argc, char* argv[]) {

	if (argc <= 1) {
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

	MatDouble a(N, N);
	MatDouble b(N, r);

	srand(time(NULL));
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			a[i][j] = rand() % 100;
		}
		b[i][0] = rand() % 40;
	}

	MatDouble L(N, N, 0.0);
	MatDouble U(N, N, 0.0);

	SystemCounterState before_state, after_state;
	clock_t begin = clock();
	try {
		before_state = getSystemCounterState();
		doolittle(a, L, U);
		for (int i = 0; i < b.n_cols(); i++) {
			VecDouble bvec(b.n_rows());
			VecDouble x(b.n_rows());
			for (int j = 0; j < b.n_rows(); j++)
				bvec[j] = b[j][i];
			solve(L, U, bvec, x);
		}
		after_state = getSystemCounterState();
	} catch (const char* e) {
		cout << e << endl;
	}
	clock_t end = clock();
	double time = (double)(end - begin);

	cout << endl << "-------------------------------------------------------------------------" << endl;
	cout << "Total number of retired instructions: " << getInstructionsRetired(before_state, after_state) << endl;
	cout << "Consumed energy by CPU in mW: " << getConsumedJoules(before_state, after_state)/(time*1000.0) << endl;
	cout << "Consumed energy by DRAM in mW: " << getDRAMConsumedJoules(before_state, after_state)/(time*1000.0) << endl;
	cout << "Bytes read from DRAM controller: " << getBytesReadFromMC(before_state, after_state) << endl;
	cout << "Bytes written to DRAM controller: " << getBytesWrittenToMC(before_state, after_state) << endl;
	cout << "Number of CPU cycles: " << getCycles(before_state, after_state) << endl;
	cout << "Instructions per cycle (IPC): " << getIPC(before_state, after_state) << endl;
	cout << "Average number of retired instructions per time interval: " << getTotalExecUsage(before_state, after_state) << endl;
	cout << "Cycles lost due to L2 cache misses but still hitting L3 cache: " << getCyclesLostDueL2CacheMisses(before_state, after_state) << endl;
	cout << "Cycles lost due to L3 cache misses: " << getCyclesLostDueL3CacheMisses(before_state, after_state) << endl;
	cout << "L2 cache misses: " << getL2CacheMisses(before_state, after_state) << endl;
	cout << "L3 cache misses: " << getL3CacheMisses(before_state, after_state) << endl;
	cout << "L2 cache hit rate: " << getL2CacheHitRatio(before_state, after_state) << endl;
	cout << "L3 cache hit rate: " << getL3CacheHitRatio(before_state, after_state) << endl;
	cout << "-------------------------------------------------------------------------" << endl;

	cout << N << " " << time * 1000.0 / CLOCKS_PER_SEC << endl;

	m->cleanup();
	return 0;
}

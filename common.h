#pragma once

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <ctime>
#include <cassert>

#include <algorithm>
#include <vector>
#include <stack>
#include <string>
#include <map>
#include <set>
using namespace std;

typedef int i32;
typedef unsigned int u32;
typedef long long i64;
typedef unsigned long long u64;

typedef double real;
typedef float f32;
typedef double f64;

const real EPS = 1e-9;

#define all(x) x.begin(), x.end()

void READ(string f) {
	freopen(f.c_str(), "r", stdin);
}

void WRITE(string f) {
	freopen(f.c_str(), "w", stdout);
}

struct STOPWATCH {
	clock_t timer;
	STOPWATCH() { timer = clock(); }
	~STOPWATCH() {
		printf("Time: %.5lfs\n", ((double)(clock() - timer)) / CLOCKS_PER_SEC);
	}
};

void debug(int n) {
    switch(n) {
        case 0: exit(0);                             // Wrong Answer
        case 1: assert(false);                       // SIGABRT
        case 2: *(int*)((u64)n - n) = 0;             // SIGSEGV
        case 3: n /= (n - 3);                        // SIGFPE
        case 4: while(1);                            // Time Limit Exceeded
        case 5: while(1) malloc(1024 * 1024);        // Memory Limit Exceeded
        case 6: malloc(32 * 1024 * 1024); while(1);  // 32MB + TLE
        case 7: while(1) puts(".");                  // Output Limit Exceeded
    }
}

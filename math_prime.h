#pragma once
#include "common.h"
#include "math_arith.h"

bool is_prime(u32 n) {
	if (n < 11) {
		return n == 2 || n == 3 || n == 5 || n == 7;
	}
	if (!(n%2) || !(n%3) || !(n%5) || !(n%7)) {
		return false;
	}
	if (n <  121) {
		return true;
	}
	for(u32 i = 11; i * i < n; i += 2) {
		if((n % i) == 0) {
			return false;
		}
	}
	return true;
}

#define GET(B,n) ( B[n >> 5] &  (1 << ( n & 31 )) )
#define SET(B,n) ( B[n >> 5] |= (1 << ( n & 31 )) )

void prime_sieve(u64 N, vector<int> &D) {
	u64 i, j, k, S = sqrt(N) + 1;
	for (i = 3; i <= S; i += 2) {
		if (!GET(D, i)) {
			for (j = i * i, k = i << 1; j < N; j += k) {
				SET(D, j);
			}
		}
	}
}

bool miller_rabin(u32 n) {
	if (n < 11) {
		return n == 2 || n == 3 || n == 5 || n == 7;
	}
	if (!(n%2) || !(n%3) || !(n%5) || !(n%7)) {
		return false;
	}
	if (n <  121) {
		return true;
	}

	u32 a[3] = {2, 7, 61};
	int i, j, s = 1;
	u32 x, d = (n - 1) >> 1;
	while(!(d & 1)) {
		d >>= 1;
		s++;
	}

	for(i = 0; i < 3; i++) {
		x = powmod(a[i], d, n);

		if(x == 1 || x == n - 1) {
			continue;
		}
		for(j = 1; j < s; j++) {
			x = ((u64)x * x) % n;
			if(x == 1) return false;
			if(x == n - 1) break;
		}
		if(j == s) {
			return false;
		}
	}
	return true;
}

bool miller_rabin(u64 n) {
	if(n < 4294967295LU) {
		return miller_rabin((u32)n);
	}

	if (!(n%2) || !(n%3) || !(n%5) || !(n%7)) {
		return false;
	}

	u64 a[10] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};
	int i, j, s = 1;
	u64 x, d = (n - 1) >> 1;
	while(!(d & 1)) {
		d >>= 1;
		s++;
	}

	for(i = 0; i < 7; i++) {
		x = powmod(a[i], d, n);

		if(x == 1 || x == n - 1) {
			continue;
		}
		for(j = 1; j < s; j++) {
			x = mulmod(x, x, n);
			if(x == 1) return false;
			if(x == n - 1) break;
		}
		if(j == s) {
			return false;
		}
	}
	return true;
}

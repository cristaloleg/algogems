#pragma once
#include "common.h"

const int MaxN = 1e5;

template<class T>
struct Segtree {
	int size;
	vector<T> D;

	Segtree(int n = 0) {
		build(n);
	}

	void build(int n) {
		size = n;
		D.resize(2 * n);
	}

	void build(int n, T* p) {
		build(n);
		for (int i = size - 1; i > 0; --i) {
			D[i] = D[i << 1] + D[i << 1 | 1];
		}
	}

	void modify(int i, T v) {
		i += size;
		D[i] = v;
		for (; i > 1; i >>= 1) {
			D[i >> 1] = D[i] + D[i ^ 1];
		}
	}

	void modify(int a, int b, T v) {
		a += size;
		b += size;
		for (; a < b; a >>= 1, b >>= 1) {
			if (a & 1) {
				D[a] += v;
				a++;
			}
			if (b & 1) {
				b--;
				D[b] += v;
			}
		}
	}

	T query(int i) {
		T res = T(0);
		for (i += size; i > 0; i >>= 1) {
			res = res + D[i];
		}
		return res;
	}

	T query(int a, int b) {
		T res = T(0);
		a += size;
		b += size;		
		for (; a < b; a >>= 1, b >>= 1) {
			if (a & 1) {
				a++;
				res = res + D[a];
			}
			if (b & 1) {
				b--;
				res = res + D[b];
			}
		}
		return res;
	}

	void push() {
		for (int i = 1; i < size; ++i) {
			D[i << 1]     = D[i << 1]     + D[i];
			D[i << 1 | 1] = D[i << 1 | 1] + D[i];
			D[i] = 0;
		}
	}
};

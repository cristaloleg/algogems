#pragma once
#include "common.h"

template<typename T>
struct Fenwick {
	int size;
	vector<T> D;

	Fenwick(int n = 0) {
		build(n);
	}

	void build(int n) {
		size = n;
		D.resize(n);
	}

	void update(int i, T v) {
		for(; i <= size; i |= i + 1) {
			D[i] += v;
		}
	}

	T query(int i) {
		T s = 0;
		for(; i >= 0; i -= i & -i) {
			s += D[i];
		}
		return s;
	}

	T query(int i, int j) {
		return query(j) - query(i - 1);
	}

	T get(int i) {
		return query(i, i);
	}

	void set(int i, T v) {
		update(i, v - get(i));
	}
};

template <class T>
struct Fenwick2D {
	int size_n, size_m;
	vector< vector<T> > D;

	Fenwick2D(int n = 0, int m = 0) {
		build(n, m);
	}

	void build(int n, int m) {
		size_n = n;
		size_m = m;
		D.resize(size_n);
		for(int i = 0; i < size_m; ++i) {
			D[i].resize(size_m);
		}
	}

	void update(int i, int j, T v = 1) {
		for(; i <= size; i |= i + 1) {
			for(int k = j; k <= size; k |= k + 1) {
				D[i][k] += v;
			}
		}
	}

	T query(int i, int j) {
		T s = 0;
		for(; i >= 0; i -= i & -i) {
			for(int k = j; k >= 0; k -= k & -k) {
				s += D[i][k];
			}
		}
		return s;
	}

	T query(int i, int j, int k, int l) { 
		return query(k, l) - query(i - 1, l) - query(k, j - 1) + query(i - 1, j - 1);
	}

	T get(int i, int j) {
		return query(i, j, i, j);
	}

	void set(int i, int j, T v) {
		update(i, j, v - get(i, j));
	}
};

template<typename T>
struct FenwickRange {
	int size;
	Fenwick<T> D, B;

	FenwickRange(int n = 0) {
		build(n);
	}

	void build(int n) {
		size = n;
		D.build(n);
		B.build(n);
	}

	void _update(int i, T v1, T v2) {
		D.update(i, v1);
		B.update(i, v2);
	}

	void update(int i, int j, T v) {
		_update(i, v, v * (1 - i));
		_update(j + 1, -v, v * b);
	}

	T _query(int i) {
		return i * D.query(i) + B.query(i);
	}

	T query(int i, int j) {
		return _query(j) - _query(i - 1);
	}

	T get(int i) {
		return query(i, i);
	}
}

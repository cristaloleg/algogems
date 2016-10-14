#pragma once
#include "common.h"

struct UnionFind {
    int size, count;
    vector<int> P, R;

    UnionFind(int n) {
    	build(n);
    }

    void build(int n) {
        size = count = n;
        P.resize(n);
        R.resize(n);
        for(int i = 0; i < n; ++i) {
        	P[i] = i;
        }
    }

    int sets() {
        return count;
    }

    int sizeofset(int i) {
        return R[i];
    }

    int find(int x) {
    	return P[x] == x ? x : P[x] = find(P[x]);
    }

    bool check(int x, int y) {
    	return find(x) == find(y);
    }

    bool unite(int x, int y) {
    	x = find(x);
    	y = find(y);
        if(x == y) {
        	return false;
        }
        if(R[x] < R[y]) { P[x] = y; R[y] += R[x]; }
        else            { P[y] = x; R[x] += R[y]; }
        count--;
        return true;
    }
};

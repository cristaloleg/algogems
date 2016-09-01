#pragma once
#include "common.h"

template<typename T>
T gcd(T a, T b) {
	return b ? gcd(a % b, b) : a;
}

template<typename T>
T lcm(T a, T b) {
	return a / gcd(a, b) * b;
}

// ax + by = gcd(a,b)
template<typename T>
T gcd_ex(T a, T b, T& x, T& y) {
	if(a == 0) {
		x = 0;
		y = 1;
		return b;
	}
	T x1, y1;
	T d = gcd_ex(b % a, a, x1, y1);
	x = y1 - (b / a) * x1;
	y = x1;
	return d;
}

// ax = 1 (mod m)
template<typename T>
T inverse(T a, T m) {
	T x, y;
	T g = gcdex(a, m, x, y);
	if (g == 1) {
		return (x % m + m) % m;
	}
	return 0;
}

template<typename T>
T mulmod(T a, T b, T m) {
    T x = 0, y = a % m;
    while(b) {
        if(b & T(1)) x = x + y;
        y = y << 1;
        if(x >= m) x = x - m;
        if(y >= m) y = y - m;
        b = b >> 1;
    }
    return x;
}

u32 powmod(u32 a, u32 n, u32 m) {
    u32 x = 1;
    while(n) {
        if(n & 1) x = ((u64)x * a) % m;
		a = ((u64)a * a) % m;
        n = n >> 1;
    }
    return x;
}

template<typename T>
T powmod(T a, T n, T m) {
    T x = 1;
    while(n) {
        if(n & T(1)) x = mulmod(x, a, m);
        a = mulmod(a, a, m);
        n = n >> 1;
    }
    return x;
}

i32 abs(i32 x) {
	i32 y = x >> 31;
	return (x + y) ^ y;
}

i64 abs(i64 x) {
	i64 y = x >> 63;
	return (x + y) ^ y;
}

i32 sign(i32 x) {
	return (x >> 31) | ((-x) >> 31);
}

i64 sign(i64 x) {
	return (x >> 63) | ((-x) >> 63);
}

i32 msb(i32 x) {
    union { double a; i32 b[2]; };
    a = x;
    return (b[1] >> 20) - 1023;
}

i32 reverse(i32 x) {
	x = ((x & 0xaaaaaaaa) >> 1) | ((x & 0x55555555) << 1);
	x = ((x & 0xcccccccc) >> 2) | ((x & 0x33333333) << 2);
	x = ((x & 0xf0f0f0f0) >> 4) | ((x & 0x0f0f0f0f) << 4);
	x = ((x & 0xff00ff00) >> 8) | ((x & 0x00ff00ff) << 8);
	x = ((x & 0xffff0000) >> 16)| ((x & 0x0000ffff) << 16);
	return x;
}

u64 X_RAND;

u64 xorshift64star() {
	X_RAND ^= X_RAND >> 12;
	X_RAND ^= X_RAND << 25;
	X_RAND ^= X_RAND >> 27;
	return X_RAND * 2685821657736338717LLU;
}

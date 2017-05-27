#pragma once
#include "common.h"

struct u128 {
	u64 hi, lo;

	u128() : lo(0), hi(0) {}
	u128(int v) : lo(v), hi(0) {}
	u128(const u128 &v) { lo = v.lo; hi = v.hi; }
	u128(const u64 &hii, const u64 &loo) { hi = hii; lo = loo; }

	int bits() const {
		int out = hi ? 64 : 0;
		u64 tmp = hi ? hi : lo;
		while(tmp) {
			tmp >>= 1;
			out++;
		}
		return out;
	}

	bool iszero() const { return !(hi || lo); }
	bool isone() const { return (hi == 0 && lo == 1); }
	operator bool() const { return (hi || lo); }
	operator int() const{ return (int)lo; }

	u128 operator +(const u128 &v) const { return u128(hi + v.hi + ((lo + v.lo) < lo), lo + v.lo); }
	u128 operator -(const u128 &v) const { return u128(hi - v.hi - ((lo - v.lo) > lo), lo - v.lo); }

	u128 operator &(const u128 &v) const { return u128(hi & v.hi, lo & v.lo); }
	u128 operator |(const u128 &v) const { return u128(hi | v.hi, lo | v.lo); }
	u128 operator ^(const u128 &v) const { return u128(hi ^ v.hi, lo ^ v.lo); }

	bool operator ==(const u128 &v) const { return (hi == v.hi) && (lo == v.lo); }
	bool operator !=(const u128 &v) const { return (hi != v.hi) || (lo != v.lo); }

	bool operator >(const u128 &v) const { return (hi == v.hi) ? (lo > v.lo) : (hi > v.hi); }
	bool operator <(const u128 &v) const { return (hi == v.hi) ? (lo < v.lo) : (hi < v.hi); }

	bool operator >=(const u128 &v) const { return (*this > v) | (*this == v); }
	bool operator <=(const u128 &v) const { return (*this < v) | (*this == v); }

	u128 operator *(const u128 &v) const {
		if(iszero() || v.isone()) {
			return *this;
		}
		if(v.iszero()) {
			return u128(0);
		}

		u128 res, tmp = v, a(*this);
		u128 zero(0), one(1);

		for (int i = 0; i < 128 && tmp; ++i) {
			if((tmp & one) != zero) {
				res = res + a;
			}
			a = a << 1;
			tmp = tmp >> 1;
		}
		return res;
	}

	/*
	inline uint128& uint128::operator*=(const uint128& b) {
		u64 a96 = hi_ >> 32;
		u64 a64 = hi_ & 0xffffffffu;
		u64 a32 = lo_ >> 32;
		u64 a00 = lo_ & 0xffffffffu;

		u64 b96 = b.hi_ >> 32;
		u64 b64 = b.hi_ & 0xffffffffu;
		u64 b32 = b.lo_ >> 32;
		u64 b00 = b.lo_ & 0xffffffffu;

		// multiply [a96 .. a00] x [b96 .. b00]
		// terms higher than c96 disappear off the high side
		// terms c96 and c64 are safe to ignore carry bit
		u64 c96 = (a96 * b00) + (a64 * b32 + a32 * b64) + (a00 * b96);
		u64 c64 = (a64 * b00) + (a32 * b32) + (a00 * b64);
		this->hi_ = (c96 << 32) + c64;
		this->lo_ = 0;

		// add terms after this one at a time to capture carry
		*this += uint128(a32 * b00) << 32;
		*this += uint128(a00 * b32) << 32;
		*this += a00 * b00;
		return *this;
	}
	*/

	void divrem(const u128 &num, const u128 &den, u128 &quo, u128 &rem) {
		if(den == u128(0)) {
			int a = 0; a = a/a;
		}
		else if(den == u128(1) || num.iszero()) {
			quo = num;
			rem = u128(0);
		}
		else if(num == den) {
			quo = u128(1);
			rem = u128(0);
		}
		else {
			quo = u128(0);
			rem = num;

			u128 copyd = den << (num.bits() - den.bits());
			u128 adder = u128(1) << (num.bits() - den.bits());

			if(copyd > rem) {
				copyd = copyd >> 1;
				adder = adder >> 1;
			}
			while(rem >= den) {
				if (rem >= copyd) {
					rem = rem - copyd;
					quo = quo | adder;
				}
				copyd = copyd >> 1;
				adder = adder >> 1;
			}
		}
	}

	u128 operator /(const u128 &v) {
		u128 rem, quo;
		divrem(*this, v, quo, rem);
		return quo;
	}

	u128 operator %(const u128 &v) {
		u128 rem, quo;
		divrem(*this, v, quo, rem);
		return rem;
	}

	u128 operator <<(const int &v) const {
		u128 res(hi, lo);
		if(v >= 128) {
			res.hi = res.lo = 0;
		} else if(v >= 64) {
			res.hi = lo << (v-64);
			res.lo = 0;
		} else {
			res.hi <<= v;
			res.hi |= lo >> (64-v);
			res.lo <<= v;
		}
		return res;
	}

	u128 operator >>(const int &v) const {
		u128 res(hi, lo);
		if(v >= 128) {
			res.hi = res.lo = 0;
		} else if(v >= 64) {
			res.lo = hi >> (v-64);
			res.hi = 0;
		} else {
			res.lo >>= v;
			res.lo |= hi << (64-v);
			res.hi >>= v;
		}
		return res;
	}

	string to_string(unsigned radix = 10) {
		if(iszero()) {
			return "0";
		}
		if(isone()) {
			return "1";
		}
		static char sz[128 + 1];
		sz[sizeof(sz) - 1] = '\0';
		u128 ii(*this), rad = u128(radix);
		int i = 128 - 1;

		while (ii != u128(0) && i) {

			u128 remainder;
			u128 copy;
			divrem(ii, rad, copy, remainder);
			sz [--i] = "0123456789"[(int)remainder];
			ii = copy;
		}
		return &sz[i];
	}

	u128(const string &S) : lo(0), hi(0) {
		if(S.empty()) {
			return;
		}
		u128 radix(10);
		string::const_iterator i = S.begin();

		if(*i == '0') {
			radix = u128(8);
			++i;
			if(i != S.end()) {
				if(*i == 'x') {
					radix = u128(16);
					++i;
				}
			}
		}

		while(i != S.end()) {
			unsigned int n;
			const char ch = *i;

			if(ch >= '0' && ch <= '9') {
				n = (ch - '0');
			} else {
				break;
			}

			(*this) = (*this) * radix;
			(*this) = (*this) + u128(n);

			++i;
		}
	}
};

template<class T, int B>
struct UExt {
	T hi, lo;

	UExt() : lo(T(0)), hi(T(0)) {}
	UExt(const T v) : hi(T(0)), lo(v) {}
	UExt(const UExt &v) { lo = v.lo; hi = v.hi; }
	UExt(const T hii, const T loo) { hi = hii; lo = loo; }

	bool iszero() const { return !(hi || lo); }
	bool isone() const { return (hi == T(0) && lo == T(1)); }
	operator bool() const { return (hi || lo); }
	operator int() const{ return (int)lo; }

	int bits() const {
		int out = hi ? B : 0;
		T tmp = hi ? hi : lo;
		while(tmp) {
			tmp = tmp >> 1;
			out++;
		}
		return out;
	}

	bool operator ==(const UExt &v) const { return (hi == v.hi) && (lo == v.lo); }
	bool operator !=(const UExt &v) const { return (hi != v.hi) || (lo != v.lo); }

	bool operator >(const UExt &v) const { return (hi == v.hi) ? (lo > v.lo) : (hi > v.hi); }
	bool operator <(const UExt &v) const { return (hi == v.hi) ? (lo < v.lo) : (hi < v.hi); }

	bool operator >=(const UExt &v) const { return (*this > v) | (*this == v); }
	bool operator <=(const UExt &v) const { return (*this < v) | (*this == v); }

	UExt operator &(const UExt &v) const { return UExt(hi & v.hi, lo & v.lo); }
	UExt operator |(const UExt &v) const { return UExt(hi | v.hi, lo | v.lo); }
	UExt operator ^(const UExt &v) const { return UExt(hi ^ v.hi, lo ^ v.lo); }

	UExt operator +(const UExt &v) const { return UExt(hi + v.hi + T((lo + v.lo) < lo), lo + v.lo); }
	UExt operator -(const UExt &v) const { return UExt(hi - v.hi - T((lo - v.lo) > lo), lo - v.lo); }

	UExt operator +(const T v) const { return UExt(hi + T((lo + v) < lo), lo + v); }
	UExt operator -(const T v) const { return UExt(hi - T((lo - v) > lo), lo - v); }

	UExt operator *(const UExt &v) const {
		if(iszero() || v.isone()) {
			return *this;
		}
		if(v.iszero()) {
			return UExt(T(0));
		}

		UExt res, tmp = v, a(*this);
		const UExt zero(T(0)), one(T(1));

		for (int i = 0; i < 2 * B && tmp; ++i) {
			if((tmp & one) != zero) {
				res = res + a;
			}
			a = a << 1;
			tmp = tmp >> 1;
		}
		return res;
	}

	static T mul(const T a , const T b, T* hi) {
		T AH = a >> (B / 2);
		T AL = (T)a;

		T BH = b >> (B / 2);
		T BL = (T)b;

		T AHBH = (T)AH * BH;
		T ALBL = (T)AL * BL;
		T AHBL = (T)AH * BL;
		T ALBH = (T)AL * BH;

		// take care of integer overflow
		T middle = AHBL + ALBH;
		if (middle < AHBL) {
			AHBH = AHBH + (1ULL << (B / 2));
		}

		T res_lo = ALBL + (middle << (B / 2));
		if (res_lo < ALBL) {
			AHBH = AHBH + (u64)1;
		}

		*hi = AHBH + middle >> (B / 2);

		return res_lo;
	}

	void divrem(const UExt &num, const UExt &den, UExt &quo, UExt &rem) {
		if(den.iszero()) {
			int a = 0; a = a/a;
		}
		else if(den.isone() || num.iszero()) {
			quo = num;
			rem = UExt(T(0));
		}
		else if(num == den) {
			quo = UExt(T(1));
			rem = UExt(T(0));
		}
		else {
			quo = UExt(T(0));
			rem = num;

			UExt copyd = den << (num.bits() - den.bits());
			UExt adder = UExt(T(1)) << (num.bits() - den.bits());

			if(copyd > rem) {
				copyd = copyd >> 1;
				adder = adder >> 1;
			}
			while(rem >= den) {
				if (rem >= copyd) {
					rem = rem - copyd;
					quo = quo | adder;
				}
				copyd = copyd >> 1;
				adder = adder >> 1;
			}
		}
	}

	UExt operator /(const UExt &v) {
		UExt rem, quo;
		divrem(*this, v, quo, rem);
		return quo;
	}

	UExt operator %(const UExt &v) {
		UExt rem, quo;
		divrem(*this, v, quo, rem);
		return rem;
	}

	UExt operator <<(const int &v) const {
		UExt res;
		if(v >= 2 * B) {
			res.hi = res.lo = 0;
		} else if(v >= B) {
			res.hi = res.lo << (v - B);
			res.lo = 0;
		} else {
			res.hi = res.hi << v;
			res.hi = res.hi | (res.lo >> (B - v));
			res.lo = res.lo << v;
		}
		return res;
	}

	UExt operator >>(const int &v) const {
		UExt res;
		if(v >= 2 * B) {
			res.hi = res.lo = 0;
		} else if(v >= B) {
			res.lo = res.hi >> (v - B);
			res.hi = 0;
		} else {
			res.lo = res.lo >> v;
			res.lo = res.lo | (res.hi << (B - v));
			res.hi = res.hi >> v;
		}
		return res;
	}

	string to_string(unsigned radix = 10) {
		if(iszero()) {
			return "0";
		}
		if(isone()) {
			return "1";
		}
		static char sz[2 * B + 1];
		sz[2 * B] = '\0';
		UExt ii(*this), rad(radix);
		int i = 2 * B - 1;

		while (!ii.iszero() && i) {

			UExt remainder;
			UExt copy;
			divrem(ii, rad, copy, remainder);
			sz [--i] = "0123456789"[(int)remainder];
			ii = copy;
		}
		return &sz[i];
	}

	UExt(const string &S) : lo(0), hi(0) {
		if(S.empty()) {
			return;
		}
		UExt radix(T(10));
		string::const_iterator i = S.begin();

		if(*i == '0') {
			radix = UExt(8);
			++i;
			if(i != S.end()) {
				if(*i == 'x') {
					radix = UExt(16);
					++i;
				}
			}
		}

		while(i != S.end()) {
			unsigned int n;
			const char ch = *i;

			if(ch >= '0' && ch <= '9') {
				n = (ch - '0');
			} else {
				break;
			}

			(*this) = (*this) * radix;
			(*this) = (*this) + UExt(n);

			++i;
		}
	}
};

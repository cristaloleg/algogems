
struct Cmplx {
	real Re, Im;
	Cmplx(real x = 0, real y = 0): Re(x), Im(y) {}
	
	real len() {
		return hypot(Re, Im);
	}

	Cmplx operator +(Cmplx &c) { return Cmplx(Re + c.Re, Im + c.Im); }
	Cmplx operator -(Cmplx &c) { return Cmplx(Re - c.Re, Im - c.Im); }

	Cmplx operator *(Cmplx &c) { 
		return Cmplx(
			Re * c.Re - Im * c.Im,
			Re * c.Im + Im * c.Re
		);
	}

	Cmplx operator /(Cmplx &c) {
		real t = c.Re * c.Re + c.Im * c.Im; 
		return Cmplx((Re * c.Re + Im * c.Im) / t,
					(Im * c.Re + Re * c.Im) / t);
	}

	// Cmplx exp(complex<ld> e) {
	// 	e = exp(e);
	// 	return Cmplx(real(e), imag(e));
	// }
};

struct Frac {
	i32 Num, Den;
	Frac(i32 a = 1, i32 b = 1): Num(a), Den(b) {}

	void normalize() {
		i32 g = gcd(Num, Den);
		Num /= g;
		Den /= g;
	}

	real toReal() {
		return (real)Num / Den;
	}

	Frac operator +(Frac &f) {
		Frac res(Num, Den);
		res += f;
		res.normalize();
		return res;
	}

	Frac operator -(Frac &f) {
		Frac res(Num, Den);
		res -= f;
		res.normalize();
		return res;
	}

	Frac operator *(Frac &f) {
		Frac res(Num, Den);
		res *= f;
		res.normalize();
		return res;
	}

	Frac operator /(Frac &f) {
		Frac res(Num, Den);
		res /= f;
		res.normalize();
		return res;
	}
};



bool is_fib(i64 n) {
	return is_square(5 * n * n + 4) || is_square(5 * n * n - 4); 
}


bool is_square(i64 x) {
	i64 goodMask = 0xC840C04048404040; // 0xC840C04048404040 computed below
	// for (int i=0; i<64; ++i) goodMask |= Long.MIN_VALUE >>> (i*i);
	
    // This tests if the 6 least significant bits are right.
    // Moving the to be tested bit to the highest position saves us masking.
    if ((goodMask << x) >= 0) {
		return false;
	}
    int numberOfTrailingZeros = 0;//Long.numberOfTrailingZeros(x);
    // Each square ends with an even number of zeros.
    if ((numberOfTrailingZeros & 1) != 0) {
		return false;
	}
    x >>= numberOfTrailingZeros;
    // Now x is either 0 or odd.
    // In binary each odd square ends with 001.
    // Postpone the sign test until now; handle zero in the branch.
    if ((x&7) != 1 | x <= 0) {
		return x == 0;
	}
    // Do it in the classical way.
    // The correctness is not trivial as the conversion from long to double is lossy!
    i64 sq = (long)sqrt(x);
    return sq * sq == x;
}

i64 isqrt(i64 x) {
	if (x <= 0) {
		return 0;
	}
	i64 s = 1, remainder = x;
	while (remainder >= s) {
		remainder -= s;
		s += 2;
	}
	return (s-1) >> 1;
}


// a^x = b (mod m)
template<typename T>
T solve(T a, T b, T m) {
	T n = (T)sqrt (m + .0) + 1;

	T an = 1;
	for (T i = 0; i < n; ++i)
		an = (an * a) % m;

	map<T,T> vals;
	for (T i = 1, cur = an; i <= n; ++i) {
		if (!vals.count(cur)) {
			vals[cur] = i;
		}
		cur = (cur * an) % m;
	}

	for (T i=0, cur=b; i<=n; ++i) {
		if (vals.count(cur)) {
			T ans = vals[cur] * n - i;
			if (ans < m) {
				return ans;
			}
		}
		cur = (cur * a) % m;
	}
	return -1;
}



#define GET2(B,n) ( B[n >> 6] &  (1 << (( n >> 1 ) & 31 )) )
#define SET2(B,n) ( B[n >> 6] |= (1 << (( n >> 1 ) & 31 )) )

#define SHIFT(k) ( 3*k*k + (1-2*k) * (1+(k&1)) )
#define STEP(k,h) ( ( (k-1) << ((k&1)^(h&1)?1:2) ) + ( (k&1) || (h&1) ? 0:2 ) )

void prime_sieve2(u32 N, vector<u32> &D) {
	u32 i, j, S = (sqrt(N) + 2) / 3 + 1;
	N = (N + 2) / 3 + 1;
	for(i = 2; i <= S; i++) {
		if(!GET2(D, i)) {
			for(j = SHIFT(i); j < N; j += STEP(i,j), j++) {
				SET2(D, j);
			}
		}
	}
}

u32 F(u32 x) {
	u32 index[] = {0,1,0,0,0};
	return ((x / 15) << 2) + index[x % 30];
}

void prime_sieve3(u32 N, vector<u32> &D) {
	u32 i, j, k;
	u32 cycle[] = {4,2,4,2,4,6,2,6};
	for(i = 0, k = 7; k < N; i++, k += cycle[i & 7]) {
		if(!GET2(D, i)) {
			for(j = F(k*k); j < N; j += k * ((j ^ 1) & 1)) {
				SET2(D, j);
			}
		}
	}
}


template<typename T>
bool miller_rabin_T(T n) {
	if(n < (T)4294967295) {
		return miller_rabin_32(n);
	} else if(n > (T)(1 << 20)) {
		return miller_rabin((u64)n);
	}

	if (!(n%T(2)) || !(n%T(3)) || !(n%T(5)) || !(n%T(7))) {
		return false;
	}

	T a[10] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};

	int i, j, r, s;
	if(n < T(4294967295LU)) {
		i = 0; r = 3;
	} else {
		i = 3; r = 10;
	}

	s = 1;
	T x, d = (n - T(1)) >> T(1);
	while(!(d & T(1))) {
		d = d >> 1;
		s++;
	}

	for(; i < r; i++) {
		x = powmod(a[i], d, n);

		if(x == T(1) || x == n - T(1)) {
			continue;
		}
		for(j = 1; j < s; j++) {
			x = (x * x) % n;
			if(x == T(1)) return false;
			if(x == n - T(1)) break;
		}
		if(j == s) {
			return false;
		}
	}
	return true;
}


i32 jacobian(i64 a, i64 n) {
	if(!a) return 0; // (0/n) = 0
	int ans = 1;
	u64 temp;
	if(a < 0) {
		a = -a;    // (a/n) = (-a/n)*(-1/n)
		if(n % 4 == 3) {
			ans = -ans; // (-1/n) = -1 if n = 3 ( mod 4 )
		} 
	}
	if(a == 1) {
		return ans; // (1/n) = 1
	}
	while(a) {
		if(a < 0) {
			a = -a;    // (a/n) = (-a/n)*(-1/n)
			if(n % 4 == 3) {
				ans = -ans;    // (-1/n) = -1 if n = 3 ( mod 4 )
			}
		}
		while(a % 2 == 0) {
			a = a / 2;    // Property (iii)
			if(n % 8 == 3 || n % 8 == 5) {
				ans = -ans;
			}
		}
		swap(a, n);    // Property (iv)
		if(a % 4 == 3 && n % 4 == 3) {
			ans = -ans; // Property (iv)
		}
		a = a % n; // because (a/p) = (a % p / p) and a % pi = (a % n) % pi if n % pi = 0
		if(a > n / 2) {
			a = a - n;
		} 
	}
	if(n == 1) {
		return ans;
	}
	return 0; 
}

i32 jacobi(i64 a, i64 m) {
	i32 ans = 1;
	a %= m;
	while(a) {
		while(!(a & 1)) {
			a >>= 1;
			if((m & 7) == 3 || (m & 7) == 5) {
				ans = -ans;
			}
		}
		swap(a, m);
		if((a & 3) == 3 && (m & 3) == 3) {
			ans = -ans;
		}
		a %= m;
	}
	if(m == 1) {
		return ans;
	}
	return 0;
}

bool solovoy_strassen(i64 n, i32 iteration) {
	if (n < 11) {
		return n == 2 || n == 3 || n == 5 || n == 7;
	}
	if (!(n%2) || !(n%3) || !(n%5) || !(n%7)) {
		return false;
	}
	if (n <  121) {
		return true;
	}

	i64 a, j, mod;
	for(i32 i = 0; i < iteration; i++) {
		a = xorshift64star() % (n - 1) + 1;
		j = (n + jacobian(a, n)) % n;
		mod = powmod(a, (n - 1) / 2, n);
		if(!j || mod != j) { 
			return false;
		}
	}
	return true;
}


template <class T, class T2>
bool lucas_selfridge(T& n, T2 unused)
{
	// сначала проверяем тривиальные случаи
	if (n < 11) {
		return n == 2 || n == 3 || n == 5 || n == 7;
	}
	if (!(n%2) || !(n%3) || !(n%5) || !(n%7)) {
		return false;
	}
	if (n <  121) {
		return true;
	}

	// проверяем, что n не является точным квадратом, иначе алгоритм даст ошибку
	// if (perfect_square(n))
		// return false;

	// алгоритм Селфриджа: находим первое число d такое, что:
	// jacobi(d,n)=-1 и оно принадлежит ряду { 5,-7,9,-11,13,... }
	T2 dd;
	for (T2 d = 5, sign = 1; ; sign = -sign, d += 2) {
		dd = d * sign;
		T g = gcd(n, d);
		if (1 < g && g < n)
			// нашли делитель - d
			return false;
		if (jacobi(T(dd), n) == -1)
			break;
	}

	// параметры Селфриджа
	T2
	p = 1,
	q = (p*p - dd) / 4;

	// разлагаем n+1 = d*2^s
	T n_1 = n;
	++n_1;
	T s, d;
	transform_num (n_1, s, d);

// 	// алгоритм Лукаса
	T
	u = 1,
	v = p,
	u2m = 1,
	v2m = p,
	qm = q,
	qm2 = q*2,
	qkd = q;

// 	for (unsigned bit = 1, bits = bits_in_number(d); bit < bits; bit++)
// 	{
// 		mulmod (u2m, v2m, n);
// 		mulmod (v2m, v2m, n);
// 		while (v2m < qm2)
// 			v2m += n;
// 		v2m -= qm2;
// 		mulmod (qm, qm, n);
// 		qm2 = qm;
// 		redouble (qm2);
// 		if (test_bit (d, bit))
// 		{
// 			T t1, t2;
// 			t1 = u2m;
// 			mulmod (t1, v, n);
// 			t2 = v2m;
// 			mulmod (t2, u, n);

// 			T t3, t4;
// 			t3 = v2m;
// 			mulmod (t3, v, n);
// 			t4 = u2m;
// 			mulmod (t4, u, n);
// 			mulmod (t4, (T)dd, n);

// 			u = t1 + t2;
// 			if (!even (u))
// 			u += n;
// 			bisect (u);
// 			u %= n;

// 			v = t3 + t4;
// 			if (!even (v))
// 				v += n;
// 			bisect (v);
// 			v %= n;
// 			mulmod (qkd, qm, n);
// 		}
// 	}

	// точно простое (или псевдо-простое)
	if (u == 0 || v == 0)
		return true;

	// довычисляем оставшиеся члены
	T qkd2 = qkd;
	redouble (qkd2);
	for (T2 r = 1; r < s; ++r)
	{
		mulmod (v, v, n);
		v -= qkd2;
		if (v < 0) v += n;
		if (v < 0) v += n;
		if (v >= n) v -= n;
		if (v >= n) v -= n;
		if (v == 0)
			return true;
		if (r < s-1)
		{
			mulmod (qkd, qkd, n);
			qkd2 = qkd;
			redouble (qkd2);
		}
	}
	return false;
}

template<typename T>
bool selfridge_pomerance(T p) {
	return ((p % 5 == 2) || (p % 5 == 3)) &&
			(powmod(T(2), p - 1, p) == 1) &&
			(fib(p + 1, p) == 0);
}


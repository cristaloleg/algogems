#pragma once
#include "common.h"
#include "math_arith.h"

real f(real x) {
	return x;
}

real simpson(real a, real b, u32 iter) {
	real x, s = 0, h = (b - a) / N;
	for (u32 i=0; i <= iter; ++i) {
		x = a + h * i;
		s += f(x) * ((i==0 || i == iter) ? 1 : ((i&1)==0) ? 2 : 4);
	}
	s *= h / 3;
	return s;
}

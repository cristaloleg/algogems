math
test gcd iter/rec speed
msb more generalized
abs, sign, reverse for i32, i64 
modular functions
fib i32, i64, T
Cmplx, Frac wrap-up
isqrt/isquare select best
prime_sieve bench for different versions
prime_sieve2,3 correctness
prime_sieve return list of primes
naive, pollard-rho factor i32, i64, T

data
unionfind iter find
segtree wrap-up
matrix
batched-list
graph structures

graph
bfs, dfs
kruskal mst


[Math]

	* math_arith.h
				gcd, lcm, binpow, powmod
				Cmplx, Frac, ContFrac
				Diophant, Pell equations
				Bit and 64-bit tricks

	* math_extint.h
				N-bit int
				prime-prob, factor(math_prime, math_factor template)

	* math_prime.h
				sieve, naive, Miller-Rabin, BPSW 

	* math_modular.h
				inverse

	* math_bigint.h
				bigint, bigfloat
				Karatsuba, FFT, Chinese Remaider Theorem(CRT)
				bigint prime-prob, bigint factor

	* math_numeric.h
				Newton, binary, ternary search
				Gauss elimination
				Diff equations
				Integration

	* math_theory.h
				Euler-phi, divisors
				Fibonacci, Catalan
				Phi-function
				modular equations
				binomial coef

	* math_factor.h
				naive, Pollard-Rho, quad-sieve

[Data]

	* data_heap.h
				sqrt, binary, binomial

	* data_fenwick.h
				Fenwick, range, 2D, 3D, kD

	* data_segtree.h
				segment tree, lazy, persistent

	* data_misc.h
				unionfind, skip list

	* data_tree.h
				Misof tree, Splay, Treap, Trie(linked, array)

	* data_matrix.h
				matrix, sparce(row-col, col-row, etc)

	* data_hash.h
				hash: map, set, multimap, multiset


GET STUFF FROM:

e-maxx:
http://e-maxx.ru/algo/
translated:
http://translate.google.com/translate?sl=auto&tl=en&js=n&prev=_t&hl=no&ie=UTF-8&u=http%3A%2F%2Fe-maxx.ru%2Falgo%2F

topcoder editorials

ARITHMETIC

issquare
issquarefree
iscube
ispower (also for bigint, find the significand as well)

fast way to calculate issquare(n), also for big n

http://stackoverflow.com/questions/295579/fastest-way-to-determine-if-an-integers-square-root-is-an-integer
http://hansliten.wordpress.com/2010/07/31/faster-square-test/
http://gmplib.org/manual/Perfect-Square-Algorithm.html#Perfect-Square-Algorithm

BIGINT

square root
n-th root
fast multiply

COMBINATORICS

generate all permutations
generate all subsets
generate all n-tuples
generate all multisets
generate all partitions
generate all set partitions

generate all subsets of submask

return next/previous (combinatorial structure of some kind)
- permutation, subset, n-tuple, partition, k-partition, k-subset, etc

evaluate binomial
evaluate multinomial

lucas' theorem

rank/unrank (combinatorial structure)
- see next/previous

COMPUTATIONAL-GEOMETRY

ccw
polygon-polygon intersection
convex hull with collinear points in O(n log n)
point in polygon that handles edge cases
point in convex polygon in O(log N)
distance from point to line, line segment, circle, polygon etc
projection of point to line (2d), plane (3d)

circle geometry http://en.wikipedia.org/wiki/Circle_segment

DATA STRUCTURE

- range minimum query in <O(n),O(log n)>, update O(log n), range update (as
  fast as possible) http://wcipeg.com/wiki/Segment_tree
- interval tree (to represent union of intervals)
- range minimum query in <O(n),O(1)> http://wcipeg.com/wiki/RMQ
- binomial heap
- fibonacci heap
- avl tree with order statistics
- splay tree
- aa tree
- scapegoat tree
- treap
- b-tree?

GAME THEORY
- retrograde analysis

GRAPH

- mincost maxflow with dijkstra (potensials)
- maxflow with dinic's algorithm
- bipartite matching with capacity >1 between each pairs such that
  the array g[][]>1 (supporting multiple of each item on each side)
- bipartite matching with edge weights=1, but where one of the sides can be
  chosen >1 times, should be faster and should use less memory than
  maxflow-tripartite
- find bridge in graph: edge that splits graph if removed
  https://secure.wikimedia.org/wikipedia/en/wiki/Bridge_%28graph_theory%29
  http://translate.googleusercontent.com/translate_c?depth=1&hl=no&ie=UTF8&prev=_t&rurl=translate.google.com&sl=auto&tl=en&u=http://e-maxx.ru/algo/bridge_searching&usg=ALkJrhgctPzHPadiea0qQgM336GtlaWm-Q
- find all k-articulation points or k-bridges: all k-tuples of (nodes, edges)
  that disconnects the graph when removed
- dilworth's theorem (minimum number of chains in a DAG or something like that)
- kth shortest path
- tarjan's scc
- feasible flow (maxflow with min/max capacity per edge)
- multicommodity flow
- out-of-kilter algorithm (but only if it isn't dominated by mincost maxflow
  with dijkstra)
- proper polynomial-time mincost maxflow
- dijkstra running in O(V log V+E) (is it of any use? consider another data
  structure than heap, for example skip-list)
- bron-kerbosch for max-clique: add vertex ordering
- online algorithms in general are strongly desired!
  - find articulation points (fast) in a graph where we add AND remove nodes
    and edges
- min-cost bipartite matching (see editorial to ipsc 2005 f "find the right
  words")

MATRIX

- matadd
- determinant
- inverse
- gaussian elimination
- gaussian elimination mod
- for square matrix A, sum A^1+A^2+...+A^n (not hard to deduce, though)
- faster algorithm for solving equation system mod, like block-wiedemann
- special algorithm for solving equation system mod 2

MISC

- chess routines (move generator, detect check, checkmate, stalemate etc)

NUMBER-THEORY

- faster powmod, without init and ugly macros
- upper bound for primtall nummer n: a(n) < n*(log n + log log n - 1/2) for n>=20
- find continued fraction from arbitrary real number
- solve pell equation
- solve more general diophantine equations
- black box factorization for huge numbers (tries a combination of appropriate
  methods such as trial division, pollard rho, elliptic curve method,
  quadratic sieve, number field sieve)
- number-theoretic transform
- cube root modulo prime http://stackoverflow.com/questions/6752374/cube-root-modulo-p-how-do-i-do-this
- chinese remainder theorem when numbers aren't pairwise coprime
  (see my solution to ipsc 2005 g "gears in action", hard input)

NUMERICAL

- laguerre's (numerical) method of finding roots of a polynomial
- robinson's recursion for solving "toeplitz" matrices
- linear programming (simplex for sparse systems)
- linear programming (interior-point method which runs in O(polynomial time))
- quadratic programming (active set)
- discrete fourier transform O(N^2)
- fast fourier transform O(N log N)
- solve equation systems
- find eigenvalues
- least squares method

PROBABILITY

given n numbers, choose k unique numbers uniformly:
http://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle
uniform random
random from gaussian distribution (box muller)
draw from other distributions

SEARCH

find kth number in array in linear time (order statistics)

STRING

- suffix array:
  for instance manber, myers O(n log n) time, O(n) memory
  kasai, lee et al O(n) lcp
  also, linear-time construction (check nong et al, 2009)
  also, O(n log^2 n) if it's short and simple (see cosmin.ro)
- aho-corasick
- given a string, calculate transition table[a][b]:
  if we have matched to pos a and see char b, transition to position in table

OTHER

take old c routines and fix them (or delete the useless ones)

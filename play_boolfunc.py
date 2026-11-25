from itertools import combinations
from collections import Counter
from math import gcd
def anf(tt):
    """ANF coefficient table from truth table (Mobius over the Boolean lattice)."""
    a = tt[:]; n = (len(tt)).bit_length() - 1
    for i in range(n):
        for m in range(len(tt)):
            if (m >> i) & 1:
                a[m] ^^= a[m ^^ (1 << i)]
    return a
ANF=anf([1,0,1,0]);print(ANF)
def hamming_weight(tt):
    return sum(tt)
def algebraic_degree(tt):
    a = anf(tt)
    return max((bin(m).count("1") for m in range(len(tt)) if a[m]), default=0)
print(algebraic_degree([1,0,1,0]))
def is_balanced(tt):
    return 2 ** hamming_weight(tt) == len(tt)
print(is_balanced([1,0,0,0]))
def walsh_spectrum(tt):
    return [sum((-1) ** (tt[x] ^ (bin(u & x).count("1") & 1))for x in range(len(tt))) for u in range(len(tt))]
print(walsh_spectrum([1,0,1]))
def nonlinearity(tt):
    return len(tt) // 2 - max(map(abs, walsh_spectrum(tt))) // 2
def annihilators(tt, max_deg):
    """All nonzero g, deg g <= max_deg, with g*tt == 0.  Brute force: n <= 5."""
    n = (len(tt)).bit_length() - 1
    masks = [m for m in range(1, 1 << n)
             if m and bin(m).count("1") <= max_deg]
    hits = []
    for r in range(1, len(masks) + 1):
        for combo in combinations(masks, r):
            g = [0] * (1 << n)
            for m in combo:
                for x in range(1 << n):
                    if bin(m & x).count("1") == bin(m).count("1"):
                        g[x] ^= 1
            if all(g[x] * tt[x] == 0 for x in range(1 << n)):
                hits.append((combo, g))
    return hits
def algebraic_immunity(tt):
    """AI(f)=min{def(g)|g in annihilators f}, via the fg=(f+1)h counting argument (AI<= ceil(n/2))."""
    neg = [1 ^ v for v in tt]
    for d in range(1, len(tt).bit_length()):
        if annihilators(tt, d) or annihilators(neg, d):
            return d
print(annihilators([1,1,1], 3));AI=algebraic_immunity([1,1,1]);print(AI<=1)

"""Exact binary matrices with columns packed into Python integers.

A matrix is a tuple of column integers. Coordinate p_i is bit i. Matrix products compose right to left.
All arithmetic uses integer XOR; there are no numerical tolerances or external packages.
"""
from itertools import permutations


def need(condition, message):
    if not condition:
        raise ValueError(message)


def identity(n):
    return tuple(1 << j for j in range(n))


def apply(M, v):
    need(type(v) is int and 0 <= v < (1 << len(M)), "Binary vector dimension")
    out = 0
    while v:
        bit = v & -v
        out ^= M[bit.bit_length() - 1]
        v ^= bit
    return out


def multiply(A, B):
    need(len(A) == len(B), "Binary matrix dimensions disagree")
    return tuple(apply(A, column) for column in B)


def product(*matrices):
    need(bool(matrices), "Empty matrix word")
    out = identity(len(matrices[0]))
    for M in matrices:
        out = multiply(out, M)
    return out


def rank(vectors):
    pivots = {}
    for v in vectors:
        while v:
            bit = v.bit_length() - 1
            if bit in pivots:
                v ^= pivots[bit]
            else:
                pivots[bit] = v
                break
    return len(pivots)


def order(M, bound):
    I = identity(len(M))
    current = I
    for k in range(1, bound + 1):
        current = multiply(M, current)
        if current == I:
            return k
    raise ValueError("Matrix order exceeds declared finite bound")


def graph_polar(n, edges):
    out = [0] * n
    need(len(set(edges)) == len(edges), "Duplicate graph edge")
    for i, j in edges:
        need(0 <= i < j < n, "Graph edge labels")
        out[i] ^= 1 << j
        out[j] ^= 1 << i
    return tuple(out)


def transvection(v, polar):
    n = len(polar)
    need(type(v) is int and 0 < v < (1 << n), "Transvection vector dimension")
    return tuple((1 << j) ^ (v if (v & polar[j]).bit_count() & 1 else 0) for j in range(n))


def polynomial_multiply(a, b):
    out = 0
    while b:
        if b & 1:
            out ^= a
        a <<= 1
        b >>= 1
    return out


def polynomial_remainder(a, b):
    need(b > 0, "Zero polynomial divisor")
    while a and a.bit_length() >= b.bit_length():
        a ^= b << (a.bit_length() - b.bit_length())
    return a


def characteristic_polynomial(M):
    """det(zI+M) over F2 using the complete Leibniz formula (n<=8)."""
    n = len(M)
    need(n <= 8, "Characteristic determinant bounded to eight dimensions")
    result = 0
    for perm in permutations(range(n)):
        term = 1
        for i, j in enumerate(perm):
            polynomial = ((M[j] >> i) & 1) | (2 if i == j else 0)
            term = polynomial_multiply(term, polynomial)
            if not term:
                break
        result ^= term
    return result


def polynomial_irreducible(poly):
    degree = poly.bit_length() - 1
    need(degree > 0, "Irreducibility needs a nonconstant polynomial")
    for d in range(1, degree // 2 + 1):
        for lower in range(1 << d):
            if polynomial_remainder(poly, (1 << d) | lower) == 0:
                return False
    return True

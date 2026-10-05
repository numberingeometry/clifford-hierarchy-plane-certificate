"""Independent dense reference for the gate U_{E6} of (C18), used by the tests.

Uses direct tensor products and Gaussian-integer matrix arithmetic, not the block formulas of
e2_hierarchy.py. Matrices are (denominator power, row-major pairs).
"""
from math import isqrt

I = (0, ((1, 0), (0, 0), (0, 0), (1, 0)))
X = (0, ((0, 0), (1, 0), (1, 0), (0, 0)))
Y = (0, ((0, 0), (0, -1), (0, 1), (0, 0)))
Z = (0, ((1, 0), (0, 0), (0, 0), (-1, 0)))
SINGLE = {"I": I, "X": X, "Y": Y, "Z": Z}
WORDS = {
    "E6": (("XII", "ZZI", "IXI", "IZZ", "IXX"), ("ZII", "ZZI")),
}


def normalize(power, values):
    values = tuple(values)
    while power > 0 and all(a % 2 == b % 2 == 0 for a, b in values):
        values = tuple((a // 2, b // 2) for a, b in values)
        power -= 1
    return power, values


def size(matrix):
    n = isqrt(len(matrix[1]))
    if n*n != len(matrix[1]):
        raise ValueError("nonsquare dense reference")
    return n


def cmul(left, right):
    a, b = left
    c, d = right
    return a*c-b*d, a*d+b*c


def identity(n):
    return 0, tuple((int(r == c), 0) for r in range(n) for c in range(n))


def product(left, right):
    n = size(left)
    if size(right) != n:
        raise ValueError("dense sizes differ")
    lhs = [[(c, left[1][r*n+c]) for c in range(n) if left[1][r*n+c] != (0, 0)]
           for r in range(n)]
    rhs = [[(c, right[1][r*n+c]) for c in range(n) if right[1][r*n+c] != (0, 0)]
           for r in range(n)]
    out = [(0, 0)] * (n*n)
    for row, entries in enumerate(lhs):
        for middle, first in entries:
            for col, second in rhs[middle]:
                a, b = cmul(first, second)
                c, d = out[row*n+col]
                out[row*n+col] = c+a, d+b
    return normalize(left[0]+right[0], out)


def adjoint(matrix):
    n = size(matrix)
    return matrix[0], tuple((matrix[1][c*n+r][0], -matrix[1][c*n+r][1])
                            for r in range(n) for c in range(n))


def tensor(left, right):
    m, n = size(left), size(right)
    return normalize(left[0]+right[0],
        (cmul(left[1][(r//n)*m+c//n], right[1][(r%n)*n+c%n])
         for r in range(m*n) for c in range(m*n)))


def pauli_word(symbols):
    out = identity(1)
    for symbol in symbols:
        out = tensor(out, SINGLE[symbol])
    return out


def pauli(label, n=5):
    if type(label) is not int or not 0 <= label < (1 << (2*n)):
        raise ValueError("reference Pauli label")
    x, z = label & ((1 << n)-1), label >> n
    # Tensor leftmost qubit is MSB, public label qubit zero is LSB.
    return pauli_word("".join("Y" if x & z & (1 << j) else
                             "X" if x & (1 << j) else
                             "Z" if z & (1 << j) else "I"
                             for j in reversed(range(n))))


def root(symbols):
    p = pauli_word(symbols)
    one = identity(size(p))
    # Independent spectral expression R=(I+P)/2 + i(I-P)/2.
    return normalize(1, ((a+c+d, a-c+d) for (a,_), (c,d) in zip(one[1], p[1])))


def word(symbols):
    out = identity(1 << len(symbols[0]))
    for symbol in symbols:
        out = product(out, root(symbol))
    return out


def archived_gate(name):
    aw, bw = WORDS[name]
    a, b = word(aw), word(bw)
    blocks = (identity(8), b, a, product(a, b))
    power = max(block[0] for block in blocks)
    out = [(0, 0)] * 1024
    for control, (exponent, data) in enumerate(blocks):
        for r in range(8):
            for c in range(8):
                real, imag = data[r*8+c]
                out[(8*control+r)*32+8*control+c] = (
                    real << (power-exponent), imag << (power-exponent))
    return normalize(power, out)


def reverse_bits(value, n):
    return sum(((value >> j) & 1) << (n-1-j) for j in range(n))


def reverse_basis(matrix):
    n = size(matrix).bit_length()-1
    dim = 1 << n
    return matrix[0], tuple(matrix[1][reverse_bits(r,n)*dim+reverse_bits(c,n)]
                            for r in range(dim) for c in range(dim))


def gate(name):
    return reverse_basis(archived_gate(name))


def conjugate(matrix, label):
    return product(product(matrix, pauli(label)), adjoint(matrix))


def support(matrix):
    n = size(matrix).bit_length()-1
    dim = 1 << n
    out = []
    for label in range(1 << (2*n)):
        x, z = label & (dim-1), label >> n
        unit = ((1,0),(0,1),(-1,0),(0,-1))[(x&z).bit_count() % 4]
        real = imag = 0
        # Hilbert-Schmidt coefficient Tr(P_label^dagger M)/dim.
        for c in range(dim):
            sign = -1 if (c&z).bit_count() % 2 else 1
            u, v = sign*unit[0], -sign*unit[1]
            a, b = cmul((u,v), matrix[1][(c^x)*dim+c])
            real += a
            imag += b
        if real or imag:
            out.append(label)
    return tuple(out)


def global_phase(matrix, phase):
    return normalize(matrix[0], (cmul(value,phase) for value in matrix[1]))


CONTROL_LABELS = (0,1,2,3,4,16,31,32,33,64,65,98,255,257,511,682,1023)


def need(value, message):
    if not value:
        raise ValueError(message)

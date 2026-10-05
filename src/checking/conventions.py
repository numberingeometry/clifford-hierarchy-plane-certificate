"""Explicit conversions between qubit-zero-low and qubit-zero-high packing."""


def reverse_qubit_bits(value, n):
    if not isinstance(value,int) or not isinstance(n,int) or n<1 or not 0<=value<1<<n:
        raise ValueError('packed coordinate outside declared dimension')
    return sum(((value>>j)&1)<<(n-1-j) for j in range(n))


def reverse_pauli_label(value,n):
    if not isinstance(value,int) or not 0<=value<1<<(2*n):
        raise ValueError('Pauli label outside declared dimension')
    return reverse_qubit_bits(value&((1<<n)-1),n)|(reverse_qubit_bits(value>>n,n)<<n)


def symplectic(a,b,n):
    mask=(1<<n)-1
    return (((a&mask)&(b>>n)).bit_count()+((b&mask)&(a>>n)).bit_count())&1


def msb_x_high_to_shared_label(value,n):
    """C16 x-high/z-low MSB coordinates <-> shared x-low/z-high LSB coordinates."""
    if not isinstance(value,int) or not 0<=value<1<<(2*n):
        raise ValueError('Pauli label outside declared dimension')
    mask=(1<<n)-1
    return reverse_qubit_bits(value>>n,n)|(reverse_qubit_bits(value&mask,n)<<n)

import typing

def intmul_2048(a, b, n):
    """Calculate a * b % n with a, b, n are 2048 bits integer"""
    assert(a < 2**2048)
    assert(b < 2**2048)
    assert(n < 2**2048)
    return (a * b) % n

def chunk(a: int) -> list[int]:
    """For 32 bits operation l = 11, chunk number every 11 bits;
       List [0] will store the LSB
       RSA2048 and RSA4096 will create poly size of 384 and 768.
    """
    xs = [None] * 384
    for i in range(384):
        xs[i] = a & 0x7FF
        a >>= 11
    return xs

def dechunk(xs: list[int]) -> int:
    a = 0
    for x in reversed(xs):
        a <<= 11
        a = a | x
    return a

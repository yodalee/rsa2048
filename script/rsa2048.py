import typing


def intmul_2048(a, b, n):
    """Calculate a * b % n with a, b, n are 2048 bits integer"""
    assert (a < 2**2048)
    assert (b < 2**2048)
    assert (n < 2**2048)
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


def CT_BFU(a, b, omega, q):
    m = (b * omega) % q
    return (a + m) % q, (a - m) % q


def GS_BFU(a, b, omega, q):
    add = (a + b) % q
    sub = ((a - b) * omega) % q
    add = add // 2 if add % 2 == 0 else (add + q) // 2
    sub = sub // 2 if sub % 2 == 0 else (sub + q) // 2
    return add, sub


def br(i: int, l: int) -> int:
    s = bin(i)[2:].zfill(l)
    return int(s[::-1], 2)


def ntt_4(xs: list[int], omega: int, prime: int) -> list[int]:
    l = len(xs)
    twiddle_factors = [pow(omega, i, prime) for i in range(l)]
    twiddle_index = 1

    distance = l//2

    while distance >= 1:
        for start in range(0, l, 2*distance):
            for i in range(distance):
                tf = twiddle_factors[br(twiddle_index, 2)]
                idx1 = start + i
                idx2 = start + i + distance
                xs[idx1], xs[idx2] = CT_BFU(xs[idx1],
                                            xs[idx2], tf, prime)
            twiddle_index += 1

        distance //= 2

    return xs

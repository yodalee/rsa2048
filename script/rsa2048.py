import typing


# Cooley–Tukey butterfly unit
def CT_BFU(a, b, omega, q):
    m = (b * omega) % q
    return (a + m) % q, (a - m) % q


# Gentleman-Sande butterfly unit
def GS_BFU(a, b, omega, q):
    add = (a + b) % q
    sub = ((a - b) * omega) % q
    add = add // 2 if add % 2 == 0 else (add + q) // 2
    sub = sub // 2 if sub % 2 == 0 else (sub + q) // 2
    return add, sub


# bit-reverse sequence generator
def br(i: int, l: int) -> int:
    s = bin(i)[2:].zfill(l)
    return int(s[::-1], 2)


# Abstract class for Rsa2048 operation
class Rsa2048:
    def __init__(self):
        self.N = 2048

    def chunk(self, a: int) -> list[int]:
        """chunk number every l bits, list [0] will store the chunk near LSB"""
        raise NotImplementedError

    def dechunk(self, l: list[int]) -> int:
        """concat list of integer of length l bits into an integer"""
        raise NotImplementedError

    def ntt(self, l: list[int]) -> list[int]:
        """Run NTT on the integer list"""
        raise NotImplementedError

    def intt(self, l: list[int]) -> list[int]:
        """Run Inverse NTT on the integer list"""
        raise NotImplementedError

    def square(self, a: int, p: int) -> int:
        """Double the montgomery form number aR mod p under modulo p
        Input: a, aR mod p
        Output: c = a^2 * R mod p

        The algorithm
        Assume that ph = NTT(chunk(p)), pm1 = NTT(chunk(p^-1 mod 2^k))
        1: ah = NTT( chunk(a) )
        2: t = dechunk( INTT(ah * ah) )
        3: th = NTT( chunk(t mod 2^k) )
        4: l = dechunk( INTT(th * pm1) )
        5: lh = NTT( chunk(l mod 2^k) )
        6: r = dechunk( INTT (lh * ph) )
        7: c = t/2k - r/2k
        8: if c < 0 then c = c + p
        9: return c
        """
        pass

    def intmul(self, a: int, b: int, p: int):
        """Multiply montgormery form number a, b, calculating a * b % p
        Input: aR mod p, bR mod p
        Output: c = abR mod p

        The algorithm
        Assume that ph = NTT(chunk(p)), pm1 = NTT(chunk(p^-1 mod 2^k))
        1: ah = NTT( chunk(a) )
        2: bh = NTT( chunk(b) )
        2: t = dechunk( INTT(ah * bh) )
        3: th = NTT( chunk(t mod 2^k) )
        4: l = dechunk( INTT(th * pm1) )
        5: lh = NTT( chunk(l mod 2^k) )
        6: r = dechunk( INTT (lh * ph) )
        7: c = t/2k - r/2k
        8: if c < 0 then c = c + p
        9: return c
        """
        pass


class Rsa2048_32b(Rsa2048):
    def chunk(self, a: int) -> list[int]:
        """32 bits chunk integer every 11 bits, into 384 chunk"""
        xs = [None] * 384
        for i in range(384):
            xs[i] = a & 0x7FF
            a >>= 11
        return xs

    def dechunk(self, xs: list[int]) -> int:
        """concat chunks of 11 bits integer"""
        a = 0
        for x in reversed(xs):
            a <<= 11
            a = a | x
        return a

    def ntt(self, l: list[int]) -> list[int]:
        assert len(l) == 384, "Invalid input to NTT method, length should be 384"
        dist = 192
        # for dist in [192, 96, 48, 24, 12, 6, 3]:
        for iter_times in range(384 // dist // 2):
            start = iter_times * dist
            zeta = 12288
            for i in range(start, start + dist):
                CT_BFU(l[i], l[i+dist], zeta, 12289)
        return l

    def intt(self, l: list[int]) -> list[int]:
        assert len(l) == 384, "Invalid input to NTT method, length should be 384"
        dist = 192
        # for dist in [192, 96, 48, 24, 12, 6, 3]:
        for iter_times in range(384 // dist // 2):
            start = iter_times * dist
            zeta = 12288
            for i in range(start, start + dist):
                GS_BFU(l[i], l[i+dist], zeta, 12289)
        return l

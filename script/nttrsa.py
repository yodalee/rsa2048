# Helper function of NTT
def CT_BFU(a, b, omega, q):
    m = (b * omega) % q
    return (a + m) % q, (a - m) % q


def GS_BFU(a, b, omega, q):
    add = (a + b) % q
    sub = ((a - b) * omega) % q
    add = add // 2 if add % 2 == 0 else (add + q) // 2
    sub = sub // 2 if sub % 2 == 0 else (sub + q) // 2
    return add, sub


# Abstract class for Rsa using NTT to speed up the multiplication
class NttRsa:
    def __init__(self, N: int, l: int, len_poly: int, q1: int, q2: int):
        assert N == 2048 or N == 4096, "N must be 2048 or 4096"
        self.N = N
        self.l = l  # number of bits in the chunk
        self.len_poly = len_poly
        self.q1 = q1
        self.q2 = q2
        self.q = q1 * q2

    def chunk(self, a: int) -> list[int]:
        """Chunk number every l bits, list [0] will store the chunk near LSB"""
        chunks = []
        mask = (1 << self.l) - 1  # Mask to extract l bits
        for _ in range(self.len_poly):
            chunks.append(a & mask)
            a >>= self.l
        return chunks

    def dechunk(self, l: list[int]) -> int:
        """Dechunk list of integers, sum up each integer by offset l bits"""
        assert len(
            l) == self.len_poly, "Input list length must match self.len_poly"
        result = 0
        for i, chunk in enumerate(l):
            result += chunk << (i * self.l)
        return result

    def crt(self, x: int, y: int) -> int:
        """Chinese Remainder Theorem
        Input: x mod q1, y mod q2
        Output: z such that z mod q1 = x, z mod q2 = y"""
        raise NotImplementedError

    def crts(self, xs: list[int], ys: list[int]) -> list[int]:
        """Apply CRT element-wise on xs and ys to produce zs
        Input: xs mod q1, ys mod q2
        Output: zs such that zs mod q1 = xs, zs mod q2 = ys
        """
        assert len(xs) == self.len_poly and len(
            ys) == self.len_poly, "Lengths of xs and ys must match len_poly"
        zs = [self.crt(x, y) for x, y in zip(xs, ys)]
        return zs

    def ntt_q1(self, l: list[int]) -> list[int]:
        """Run NTT on the integer list"""
        raise NotImplementedError

    def intt_q1(self, l: list[int]) -> list[int]:
        """Run Inverse NTT on the integer list"""
        raise NotImplementedError

    def ntt_q2(self, l: list[int]) -> list[int]:
        """Run NTT on the integer list"""
        raise NotImplementedError

    def intt_q2(self, l: list[int]) -> list[int]:
        """Run Inverse NTT on the integer list"""
        raise NotImplementedError

    def mul_q1(self, a: list[int], b: list[int]) -> list[int]:
        """Multiply two NTT form numbers a, b under modulo q1
        Input: a, b in NTT form
        Output: c = ab in NTT form
        """
        raise NotImplementedError

    def mul_q2(self, a: list[int], b: list[int]) -> list[int]:
        """Multiply two NTT form numbers a, b under modulo q2
        Input: a, b in NTT form
        Output: c = ab in NTT form
        """
        raise NotImplementedError

    def square(self, a: int,
               ph1: list[int], ph2: list[int],
               pm1: list[int], pm2: list[int]) -> int:
        """Square the montgomery form number aR mod p under modulo p
        Input: a, aR mod p
        Output: c = a^2 * R mod p

        The algorithm
        Assume that ph1,ph2 = NTT(chunk(p)), pm1,pm2 = NTT(chunk(p^-1 mod 2^k))
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

        # Calculate t = a * a
        al = self.chunk(a)
        ah1 = self.ntt_q1(al)
        ah2 = self.ntt_q2(al)
        sqrh1 = self.mul_q1(ah1, ah1)
        sqrh2 = self.mul_q2(ah2, ah2)
        sqr1 = self.intt_q1(sqrh1)
        sqr2 = self.intt_q2(sqrh2)
        c = self.crts(sqr1, sqr2)
        t = self.dechunk(c)

        # l = (t mod 2**2048) * minpinv
        c_low = self.lower(c)
        ch1 = self.ntt_q1(c_low)
        ch2 = self.ntt_q2(c_low)
        lh1 = self.mul_q1(ch1, pm1)
        lh2 = self.mul_q2(ch2, pm2)
        l1 = self.intt_q1(lh1)
        l2 = self.intt_q2(lh2)
        l = self.crts(l1, l2)

        # lp = l * P
        l_low = self.lower(l)
        lh1 = self.ntt_q1(l_low)
        lh2 = self.ntt_q2(l_low)
        lph1 = self.mul_q1(lh1, ph1)
        lph2 = self.mul_q2(lh2, ph2)
        lp1 = self.intt_q1(lph1)
        lp2 = self.intt_q2(lph2)
        lp = self.crts(lp1, lp2)
        lp = self.dechunk(lp)

        return self.addhigh(t, lp)

    def intmul(self, a: int, b: int,
               ph1: list[int], ph2: list[int],
               pm1: list[int], pm2: list[int]) -> int:
        """Multiply montgormery form number a, b, calculating a * b % p
        Input: aR mod p, bR mod p
        Output: c = abR mod p

        The algorithm
        Assume that ph1, ph2 = NTT(chunk(p)), pm1, pm2 = NTT(chunk(p^-1 mod 2^k))
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
        # t = a * b
        al = self.chunk(a)
        bl = self.chunk(b)
        ah1 = self.ntt_q1(al)
        ah2 = self.ntt_q2(al)
        bh1 = self.ntt_q1(bl)
        bh2 = self.ntt_q2(bl)
        abh1 = self.mul_q1(ah1, bh1)
        abh2 = self.mul_q2(ah2, bh2)
        ab1 = self.intt_q1(abh1)
        ab2 = self.intt_q2(abh2)
        ab = self.crts(ab1, ab2)
        t = self.dechunk(ab)

        # l = (t mod 2**2048) * minpinv
        c_low = self.lower(ab)
        ch1 = self.ntt_q1(c_low)
        ch2 = self.ntt_q2(c_low)
        lh1 = self.mul_q1(ch1, pm1)
        lh2 = self.mul_q2(ch2, pm2)
        l1 = self.intt_q1(lh1)
        l2 = self.intt_q2(lh2)
        l = self.crts(l1, l2)

        # lp = l * P
        l_low = self.lower(l)
        lh1 = self.ntt_q1(l_low)
        lh2 = self.ntt_q2(l_low)
        lph1 = self.mul_q1(lh1, ph1)
        lph2 = self.mul_q2(lh2, ph2)
        lp1 = self.intt_q1(lph1)
        lp2 = self.intt_q2(lph2)
        lp = self.crts(lp1, lp2)
        lp = self.dechunk(lp)

        return self.addhigh(t, lp)
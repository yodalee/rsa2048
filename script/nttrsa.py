import math

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
        # Modulus p and their derivative, initially None
        self.p = None
        self.ph1 = None
        self.ph2 = None
        self.pm1 = None
        self.pm2 = None
        # To convert to montgomery form
        self.r = None
        self.rsqr = None

    def setp(self, p: int):
        """Set the modulus p and calculate derived fields."""
        self.p = p
        # derive NTT form of p
        p_chunks = self.chunk(p)
        self.ph1 = self.ntt_q1(p_chunks)  # NTT(chunk(p)) for q1
        self.ph2 = self.ntt_q2(p_chunks)  # NTT(chunk(p)) for q2
        # derive NTT form of p^-1 mod 2^N
        p_inv = self.qinv(p)
        p_inv = pow(p, -1, 1 << 2048)  # p^-1 mod q
        p_inv_chunks = self.chunk(p_inv)
        self.pm1 = self.ntt_q1(p_inv_chunks)  # NTT(chunk(p^-1 mod 2^k)) for q1
        self.pm2 = self.ntt_q2(p_inv_chunks)  # NTT(chunk(p^-1 mod 2^k)) for q2
        # derive R^2 mod p for conversion to montgomery form
        self.r = (1 << self.N) % p
        self.rsqr = pow(self.r, 2, p)  # R^2 mod p

    def qinv(self, p: int) -> int:
        """Calculate p^-1 mod 2^N by hensel_lifting"""
        # initial setting x*y - 1 = 2x-x^2-1 = -(x-1)^2
        # since p is odd, (x-1)^2 = 4*k
        r = 1 << self.N
        x = p & 0xFFFFFFFF
        y = (2 - x) % r

        for _ in range(int(math.log2(self.N // 2))):
            y = (2 * y - p * y * y) % r
        return y

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

    def lower(self, l: list[int]) -> list[int]:
        """Extract the lower part of a chunked number"""
        assert len(
            l) == self.len_poly, "Input list length must match self.len_poly"
        ret = [0] * self.len_poly
        mask = (1 << self.l) - 1
        remain = 0
        n_chunks = (self.N + self.l - 1) // self.l

        for i, chunk in enumerate(l[:n_chunks]):
            sum = (chunk + remain)
            ret[i] = sum & mask
            remain = sum >> self.l
        # Modulo final block to N bits
        ret[n_chunks-1] &= ((1 << self.N % self.l) - 1)
        return ret

    def square(self, a: int) -> int:
        """Square the montgomery form number aR mod p under modulo p
        Input: a, aR mod p
        Output: c = a^2 * R mod p

        It follows the montgomery reduction algorithm that:
          Mont(a, b) = high(a * b + p * low(p^-1 * low(a * b)))
        The R used here is 2^self.N

        The algorithm
        Assume that ph1,ph2 = NTT(chunk(p)), pm1,pm2 = NTT(chunk(p^-1 mod 2^k))
        1: ah = NTT( chunk(a) )
        2: t = dechunk( INTT(ah * ah) )
        3: th = NTT( chunk(t mod R) )
        4: l = dechunk( INTT(th * pm1) )
        5: lh = NTT( chunk(l mod R) )
        6: r = dechunk( INTT (lh * ph) )
        7: c = t/R - r/R
        8: if c < 0 then c = c + p
        9: return c
        """
        if self.p is None:
            raise ValueError("Modulus p has not been set. Call setp() first.")

        # Calculate t = a * a
        al = self.chunk(a)
        ah1 = self.ntt_q1(al)
        ah2 = self.ntt_q2(al)
        sqrh1 = self.mul_q1(ah1, ah1)
        sqrh2 = self.mul_q2(ah2, ah2)
        sqrl1 = self.intt_q1(sqrh1)
        sqrl2 = self.intt_q2(sqrh2)
        sqrl = self.crts(sqrl1, sqrl2)
        t = self.dechunk(sqrl)

        # l = (t mod R) * minpinv
        t_lowl = self.lower(sqrl)
        th1 = self.ntt_q1(t_lowl)
        th2 = self.ntt_q2(t_lowl)
        lh1 = self.mul_q1(th1, self.pm1)
        lh2 = self.mul_q2(th2, self.pm2)
        l1 = self.intt_q1(lh1)
        l2 = self.intt_q2(lh2)
        ll = self.crts(l1, l2)

        # lp = l * p
        l_lowl = self.lower(ll)
        lh1 = self.ntt_q1(l_lowl)
        lh2 = self.ntt_q2(l_lowl)
        lph1 = self.mul_q1(lh1, self.ph1)
        lph2 = self.mul_q2(lh2, self.ph2)
        lpl1 = self.intt_q1(lph1)
        lpl2 = self.intt_q2(lph2)
        lpl = self.crts(lpl1, lpl2)
        lp = self.dechunk(lpl)

        high = (t >> self.N) - (lp >> self.N)
        if high < 0:
            high += self.p
        return high

    def multiply(self, a: int, b: int) -> int:
        """Multiply montgormery form number a, b, calculating a * b % p
        Input: aR mod p, bR mod p
        Output: c = abR mod p

        It follows the montgomery reduction algorithm that:
          Mont(a, b) = high(a * b + p * low(p^-1 * low(a * b)))
        The R used here is 2^self.N

        The algorithm
        Assume that ph1, ph2 = NTT(chunk(p)), pm1, pm2 = NTT(chunk(p^-1 mod 2^k))
        1: ah = NTT( chunk(a) )
        2: bh = NTT( chunk(b) )
        2: t = dechunk( INTT(ah * bh) )
        3: th = NTT( chunk(t mod R) )
        4: l = dechunk( INTT(th * pm) )
        5: lh = NTT( chunk(l mod R) )
        6: r = dechunk( INTT (lh * ph) )
        7: c = t/R - r/R
        8: if c < 0 then c = c + p
        9: return c
        """
        if self.p is None:
            raise ValueError("Modulus p has not been set. Call setp() first.")

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
        abl = self.crts(ab1, ab2)
        t = self.dechunk(abl)

        # l = (t mod R) * minpinv
        t_lowl = self.lower(abl)
        th1 = self.ntt_q1(t_lowl)
        th2 = self.ntt_q2(t_lowl)
        lh1 = self.mul_q1(th1, self.pm1)
        lh2 = self.mul_q2(th2, self.pm2)
        l1 = self.intt_q1(lh1)
        l2 = self.intt_q2(lh2)
        ll = self.crts(l1, l2)

        # lp = l * P
        l_lowl = self.lower(ll)
        lh1 = self.ntt_q1(l_lowl)
        lh2 = self.ntt_q2(l_lowl)
        lph1 = self.mul_q1(lh1, self.ph1)
        lph2 = self.mul_q2(lh2, self.ph2)
        lp1 = self.intt_q1(lph1)
        lp2 = self.intt_q2(lph2)
        lpl = self.crts(lp1, lp2)
        lp = self.dechunk(lpl)

        # c = t - lp
        high = (t >> self.N) - (lp >> self.N)
        if high < 0:
            high += self.p
        return high

    def expmod_public(self, a: int, e: int) -> int:
        """Exponentiate a to the power of e under modulo p
        Input: a, e
        Output: c = a^e mod p

        The algorithm
        1. Convert a to montgomery form
        2. Square and multiply
        3. Convert back to normal form
        """
        if self.p is None:
            raise ValueError("Modulus p has not been set. Call setp() first.")
        if e < 0 or e >= (1 << 32):
            raise ValueError("Exponent e must be positive and less than 2^32")
        if e == 0:
            return 1

        # Convert a to montgomery form
        monta = self.multiply(a, self.rsqr)

        c = self.r
        binary = [int(d) for d in bin(e)[2:]]
        for b in binary:
            c = self.square(c)
            if b == 1:
                c = self.multiply(c, monta)

        # Convert back to normal form
        c = self.multiply(c, 1)
        return c

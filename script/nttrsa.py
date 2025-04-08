
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
        """chunk number every l bits, list [0] will store the chunk near LSB"""
        raise NotImplementedError

    def dechunk(self, l: list[int]) -> int:
        """dechunk list of integer, sum up each integer by offset l bits"""
        raise NotImplementedError

    def crt(self, a: int, b: int) -> int:
        """Chinese Remainder Theorem
        Input: a mod q1, b mod q2
        Output: x such that x mod q1 = a, x mod q2 = b"""
        raise NotImplementedError

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

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
    def __init__(self):
        Rsa2048.__init__(self)
        # 1 is the placeholder
        self.zetas_12289 = [
            1, 81, 6561, 3014, 10643, 1853, 2625, 3712,
            5736, 9923, 4978, 9970, 8785, 11112, 2975, 7484,
            4043, 7969, 6461, 7203, 5860, 7678, 7468, 2747,
            1305, 7393, 8961, 790, 2545, 9521, 9283, 2294,
            1479, 9198, 7698, 9088, 11077, 140, 11340, 9154,
            4134, 3051, 1351, 11119, 3542, 4255, 563, 8736,
            7143, 1000, 7266, 10963, 3195, 726, 9650, 7443,
            722, 9326, 5777, 955, 3621, 10654, 2744, 1062,
            12288, 12208, 5728, 9275, 1646, 10436, 9664, 8577,
            6553, 2366, 7311, 2319, 3504, 1177, 9314, 4805,
            8246, 4320, 5828, 5086, 6429, 4611, 4821, 9542,
            10984, 4896, 3328, 11499, 9744, 2768, 3006, 9995,
            10810, 3091, 4591, 3201, 1212, 12149, 949, 3135,
            8155, 9238, 10938, 1170, 8747, 8034, 11726, 3553,
            5146, 11289, 5023, 1326, 9094, 11563, 2639, 4846,
            11567, 2963, 6512, 11334, 8668, 1635, 9545, 11227]

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
        zeta_idx = 1
        for dist in [192, 96, 48, 24, 12, 6, 3]:
            for iter_times in range(384 // dist // 2):
                start = iter_times * dist * 2
                idx = br(zeta_idx, 7)
                zeta = self.zetas_12289[idx]
                for i in range(start, start + dist):
                    l[i], l[i+dist] = CT_BFU(l[i], l[i+dist], zeta, 12289)
                zeta_idx += 1
        return l

    def intt(self, l: list[int]) -> list[int]:
        assert len(l) == 384, "Invalid input to NTT method, length should be 384"
        zeta_idx = 127
        for dist in [3, 6, 12, 24, 48, 96, 192]:  # , 96, 48, 24, 12, 6, 3]:
            for iter_times in range(384 // dist // 2):
                start = iter_times * dist * 2
                idx = br(zeta_idx, 7)
                zeta = self.zetas_12289[idx]
                for i in range(start, start + dist):
                    l[i], l[i+dist] = GS_BFU(l[i], l[i+dist], zeta, 12289)
                zeta_idx -= 1
        return l

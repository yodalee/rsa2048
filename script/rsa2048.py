from nttrsa import NttRsa, CT_BFU, GS_BFU


class NttRsa2048_32b(NttRsa):
    """
    NTT-RSA 2048-bit key size with on 32-bit processor
    """

    def __init__(self):
        super().__init__(2048, 11, 384, 12289, 65537)
        self.ntt_len = 128
        self.ntt_round = 7
        self.ntt_index = [
            1, 33, 17, 49, 9,  41, 25, 57, 5, 37, 21, 53, 13, 45, 29, 61,
            3, 35, 19, 51, 11, 43, 27, 59, 7, 39, 23, 55, 15, 47, 31, 63,
            2, 34, 18, 50, 10, 42, 26, 58, 6, 38, 22, 54, 14, 46, 30, 62,
            4, 36, 20, 52, 12, 44, 28, 60, 8, 40, 24, 56, 16, 48, 32, 64
        ]

        self.zetas1 = [
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
            11567, 2963, 6512, 11334, 8668, 1635, 9545, 11227
        ]
        self.zetas2 = [
            1, 4938, 4080, 27181, 2, 9876, 8160, 54362,
            4, 19752, 16320, 43187, 8, 39504, 32640, 20837,
            16, 13471, 65280, 41674, 32, 26942, 65023, 17811,
            64, 53884, 64509, 35622, 128, 42231, 63481, 5707,
            256, 18925, 61425, 11414, 512, 37850, 57313, 22828,
            1024, 10163, 49089, 45656, 2048, 20326, 32641, 25775,
            4096, 40652, 65282, 51550, 8192, 15767, 65027, 37563,
            16384, 31534, 64517, 9589, 32768, 63068, 63497, 19178,
            65536, 60599, 61457, 38356, 65535, 55661, 57377, 11175,
            65533, 45785, 49217, 22350, 65529, 26033, 32897, 44700,
            65521, 52066, 257, 23863, 65505, 38595, 514, 47726,
            65473, 11653, 1028, 29915, 65409, 23306, 2056, 59830,
            65281, 46612, 4112, 54123, 65025, 27687, 8224, 42709,
            64513, 55374, 16448, 19881, 63489, 45211, 32896, 39762,
            61441, 24885, 255, 13987, 57345, 49770, 510, 27974,
            49153, 34003, 1020, 55948, 32769, 2469, 2040, 46359,
        ]

    def crt(self, x: int, y: int) -> int:
        """Chinese Remainder Theorem
        Input: x mod q1, y mod q2
        Output: z such that z mod q1 = x, z mod q2 = y"""
        # qinv = 45373 = 12289 ** 65535 % 65537
        q1inv = 45373
        return x + (((y-x) * q1inv % self.q2) * self.q1) % self.q

    def ntt(self, l: list[int], zetas: list[int], q: int) -> list[int]:
        assert len(
            l) == self.len_poly, f"NTT: Length of input list must be {self.len_poly}"
        assert len(
            zetas) == self.ntt_len, f"NTT: Length of zetas must be {self.ntt_len}"

        for i, dist in enumerate([192, 96, 48, 24, 12, 6, 3]):
            zeta_idx = -1 * (1 << i)
            for start in range(0, len(l), dist * 2):
                zeta = zetas[self.ntt_index[zeta_idx]]
                zeta_idx += 1

                for j in range(dist):
                    # print(f"NTT: {start + j} {start + j + dist} {zeta}")
                    l[start + j], l[start + j + dist] = CT_BFU(
                        l[start + j], l[start + j + dist], zeta, q)
        return l

    def intt(self, l: list[int], zetas: list[int], q: int) -> list[int]:
        assert len(
            l) == self.len_poly, f"intt: Length of input list must be {self.len_poly}"
        assert len(
            zetas) == self.ntt_len, f"intt: Length of zetas must be {self.ntt_len}"

        for i, dist in enumerate([3, 6, 12, 24, 48, 96, 192]):
            zeta_idx = -1 * (1 << (6-i))
            for start in range(0, len(l), dist * 2):
                idx = self.ntt_index[zeta_idx]
                zeta = zetas[128-idx]
                zeta_idx += 1

                for j in range(dist):
                    # print(f"INTT: {start + j} {start + j + dist} {zeta}")
                    l[start + j], l[start + j + dist] = GS_BFU(
                        l[start + j], l[start + j + dist], zeta, q)
        return l

    def ntt_q1(self, l: list[int]) -> list[int]:
        """Run NTT on the integer list"""
        return self.ntt(l[:], self.zetas1, self.q1)

    def intt_q1(self, l: list[int]) -> list[int]:
        """Run Inverse NTT on the integer list"""
        return self.intt(l[:], self.zetas1, self.q1)

    def ntt_q2(self, l: list[int]) -> list[int]:
        """Run NTT on the integer list"""
        return self.ntt(l[:], self.zetas2, self.q2)

    def intt_q2(self, l: list[int]) -> list[int]:
        """Run Inverse NTT on the integer list"""
        return self.intt(l[:], self.zetas2, self.q2)

    def mul_q1(self, a: list[int], b: list[int]) -> list[int]:
        assert len(
            a) == self.len_poly, f"mul_q1: Length of input list a must be {self.len_poly}"
        assert len(
            b) == self.len_poly, f"mul_q1: Length of input list b must be {self.len_poly}"
        c = [0] * self.len_poly

        # Multiply a2 x^2 + a1 x + a0 with b2 x^2 + b1 x + b0 Under NTT domain of x^3 - omega
        for i in range(self.ntt_len):
            idx = self.ntt_index[i//2] + (64 if i % 2 == 1 else 0)
            omega = self.zetas1[idx % 128]
            a0 = a[3*i + 0]
            a1 = a[3*i + 1]
            a2 = a[3*i + 2]
            b0 = b[3*i + 0]
            b1 = b[3*i + 1]
            b2 = b[3*i + 2]
            # c0 = a0b0 + omega(a2b1 +a1b2)
            # c1 = a1b0 + a0b1 + omega(a2b2)
            # c2 = a2b0 + a1b1 + a0b2
            c[3*i + 0] = (a0 * b0 + omega * (a2 * b1 + a1 * b2)) % self.q1
            c[3*i + 1] = (a1 * b0 + a0 * b1 + omega * (a2 * b2)) % self.q1
            c[3*i + 2] = (a2 * b0 + a1 * b1 + a0 * b2) % self.q1
        return c

    def mul_q2(self, a: list[int], b: list[int]) -> list[int]:
        assert len(
            a) == self.len_poly, f"mul_q2: Length of input list a must be {self.len_poly}"
        assert len(
            b) == self.len_poly, f"mul_q2: Length of input list b must be {self.len_poly}"

        c = [0] * self.len_poly

        # Multiply a2 x^2 + a1 x + a0 with b2 x^2 + b1 x + b0 Under NTT domain of x^3 - omega
        for i in range(self.ntt_len):
            idx = self.ntt_index[i//2] + (64 if i % 2 == 1 else 0)
            omega = self.zetas2[idx % 128]
            a0 = a[3*i + 0]
            a1 = a[3*i + 1]
            a2 = a[3*i + 2]
            b0 = b[3*i + 0]
            b1 = b[3*i + 1]
            b2 = b[3*i + 2]
            # c0 = a0b0 + omega(a2b1 +a1b2)
            # c1 = a1b0 + a0b1 + omega(a2b2)
            # c2 = a2b0 + a1b1 + a0b2
            c[3*i + 0] = (a0 * b0 + omega * (a2 * b1 + a1 * b2)) % self.q2
            c[3*i + 1] = (a1 * b0 + a0 * b1 + omega * (a2 * b2)) % self.q2
            c[3*i + 2] = (a2 * b0 + a1 * b1 + a0 * b2) % self.q2
        return c

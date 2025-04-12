import unittest
from rsa2048 import *


class TestZetasGen(unittest.TestCase):
    def test_zetas_gen_output(self):
        gen = ZetasGen()
        zetas = []
        for _ in range(127):
            zetas.append(next(gen))
        self.assertEqual(zetas,
                         [
                             # R0
                             64,
                             # R1
                             32, 64,
                             # R2
                             16, 48, 32, 64,
                             # R3
                             8, 40, 24, 56, 16, 48, 32, 64,
                             # R4
                             4, 36, 20, 52, 12, 44, 28, 60, 8, 40, 24, 56, 16, 48, 32, 64,
                             # R5
                             2, 34, 18, 50, 10, 42, 26, 58, 6, 38, 22, 54, 14, 46, 30, 62,
                             4, 36, 20, 52, 12, 44, 28, 60, 8, 40, 24, 56, 16, 48, 32, 64,
                             # R6
                             1, 33, 17, 49, 9,  41, 25, 57, 5, 37, 21, 53, 13, 45, 29, 61,
                             3, 35, 19, 51, 11, 43, 27, 59, 7, 39, 23, 55, 15, 47, 31, 63,
                             2, 34, 18, 50, 10, 42, 26, 58, 6, 38, 22, 54, 14, 46, 30, 62,
                             4, 36, 20, 52, 12, 44, 28, 60, 8, 40, 24, 56, 16, 48, 32, 64,
                         ])
        self.assertEqual(len(zetas), 127)
        # Add further checks if specific properties of the output list need to be validated


class TestNttRsa2048_32b(unittest.TestCase):
    def test_ntt_const(self):
        rsa = NttRsa2048_32b()
        xs = [3 if i == 0 else 0 for i in range(384)]
        ys = rsa.ntt_q1(xs)
        for i in range(384):
            self.assertEqual(ys[i], 3 if i % 3 == 0 else 0)

    def test_ntt_x3(self):
        rsa = NttRsa2048_32b()
        xs = [1 if i == 3 else 0 for i in range(384)]
        ys = rsa.ntt_q1(xs)
        for i, p in enumerate([1, 33, 17, 49, 9,  41, 25, 57, 5, 37, 21, 53, 13, 45, 29, 61,
                               3, 35, 19, 51, 11, 43, 27, 59, 7, 39, 23, 55, 15, 47, 31, 63,
                               2, 34, 18, 50, 10, 42, 26, 58, 6, 38, 22, 54, 14, 46, 30, 62,
                               4, 36, 20, 52, 12, 44, 28, 60, 8, 40, 24, 56, 16, 48, 32, 64]):
            self.assertEqual(ys[i*6], pow(81, p, 12289))
            self.assertEqual(ys[i*6+1], 0)
            self.assertEqual(ys[i*6+2], 0)
            self.assertEqual(ys[i*6+3], pow(81, p+64, 12289))
            self.assertEqual(ys[i*6+4], 0)
            self.assertEqual(ys[i*6+5], 0)


if __name__ == '__main__':
    unittest.main()

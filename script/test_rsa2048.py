import unittest
from rsa2048 import *


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
        for i, p in enumerate(rsa.ntt_index):
            self.assertEqual(ys[i*6], pow(81, p, 12289))
            self.assertEqual(ys[i*6+1], 0)
            self.assertEqual(ys[i*6+2], 0)
            self.assertEqual(ys[i*6+3], pow(81, p+64, 12289))
            self.assertEqual(ys[i*6+4], 0)
            self.assertEqual(ys[i*6+5], 0)

    def test_intt_const(self):
        rsa = NttRsa2048_32b()
        xs = [3 if i % 3 == 0 else 0 for i in range(384)]
        ys = rsa.intt_q1(xs)
        for i in range(384):
            self.assertEqual(ys[i], 3 if i == 0 else 0)

    def test_intt_x3(self):
        rsa = NttRsa2048_32b()
        xs = [0] * 384
        for i, p in enumerate(rsa.ntt_index):
            xs[i*6+0] = pow(81, p, 12289)
            xs[i*6+3] = pow(81, p+64, 12289)
        ys = rsa.intt_q1(xs)
        for i in range(384):
            self.assertEqual(ys[i], 1 if i == 3 else 0)


if __name__ == '__main__':
    unittest.main()

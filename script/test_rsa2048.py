import unittest
from rsa2048 import *
import random


class TestNttRsa2048_32b(unittest.TestCase):
    def setUp(self):
        self.rsa = NttRsa2048_32b()

    def test_ntt_const(self):
        xs = [3 if i == 0 else 0 for i in range(384)]
        ys = self.rsa.ntt_q1(xs)
        for i in range(384):
            self.assertEqual(ys[i], 3 if i % 3 == 0 else 0)

    def test_ntt_x3(self):
        xs = [1 if i == 3 else 0 for i in range(384)]
        ys = self.rsa.ntt_q1(xs)
        for i, p in enumerate(self.rsa.ntt_index):
            self.assertEqual(ys[i*6], pow(81, p, 12289))
            self.assertEqual(ys[i*6+1], 0)
            self.assertEqual(ys[i*6+2], 0)
            self.assertEqual(ys[i*6+3], pow(81, p+64, 12289))
            self.assertEqual(ys[i*6+4], 0)
            self.assertEqual(ys[i*6+5], 0)

    def test_intt_const(self):
        xs = [3 if i % 3 == 0 else 0 for i in range(384)]
        ys = self.rsa.intt_q1(xs)
        for i in range(384):
            self.assertEqual(ys[i], 3 if i == 0 else 0)

    def test_intt_x3(self):
        xs = [0] * 384
        for i, p in enumerate(self.rsa.ntt_index):
            xs[i*6+0] = pow(81, p, 12289)
            xs[i*6+3] = pow(81, p+64, 12289)
        ys = self.rsa.intt_q1(xs)
        for i in range(384):
            self.assertEqual(ys[i], 1 if i == 3 else 0)

    def test_multiply_q1(self):
        # a = 1 + 2x + 3x^2 + ... + 384x^383
        a = [(i + 1) for i in range(384)]
        # b = 384 + 383x + 382x^2 + ... + 1x^383
        b = [(384 - i) for i in range(384)]

        # Schoolbook multiplication
        c_schoolbook = [0] * 384
        for i in range(384):
            for j in range(384):
                if i + j < 384:
                    c_schoolbook[i + j] += a[i] * b[j]
                else:
                    c_schoolbook[i + j - 384] += a[i] * b[j]

        # NTT-based multiplication
        a_ntt = self.rsa.ntt_q1(a)
        b_ntt = self.rsa.ntt_q1(b)
        c_ntt = self.rsa.mul_q1(a_ntt, b_ntt)
        c_ntt_result = self.rsa.intt_q1(c_ntt)

        # Check if both methods produce the same result
        for i in range(384):
            self.assertEqual(c_schoolbook[i] % self.rsa.q1, c_ntt_result[i])

    def test_multiply_q2(self):
        # a = 1 + 2x + 3x^2 + ... + 384x^383
        a = [(i + 1) for i in range(384)]
        # b = 384 + 383x + 382x^2 + ... + 1x^383
        b = [(384 - i) for i in range(384)]

        # Schoolbook multiplication
        c_schoolbook = [0] * 384
        for i in range(384):
            for j in range(384):
                if i + j < 384:
                    c_schoolbook[i + j] += a[i] * b[j]
                else:
                    c_schoolbook[i + j - 384] += a[i] * b[j]

        # NTT-based multiplication
        a_ntt = self.rsa.ntt_q2(a)
        b_ntt = self.rsa.ntt_q2(b)
        c_ntt = self.rsa.mul_q2(a_ntt, b_ntt)
        c_ntt_result = self.rsa.intt_q2(c_ntt)

        # Check if both methods produce the same result
        for i in range(384):
            self.assertEqual(c_schoolbook[i] % self.rsa.q2, c_ntt_result[i])


if __name__ == '__main__':
    unittest.main()

import unittest
import random
from rsa2048 import *


class TestRsa2048_32b(unittest.TestCase):
    def setUp(self):
        self.rsa = Rsa2048_32b()

    @unittest.skip
    def test_intmul_2048(self):
        # Generate random 2048-bit integers a, b, and n
        a = random.getrandbits(2048)
        b = random.getrandbits(2048)
        n = random.getrandbits(2048)

        # Ensure n is odd (for RSA encryption)
        if n % 2 == 0:
            n += 1

        # Calculate the result using the function
        result = self.rsa.intmul(a, b, n)

        # Check if the result is equal to 1
        self.assertEqual(result, (a * b) %
                         n, "Result of intmul_2048 is incorrect")

    def test_chunk_dechunk(self):
        a = random.getrandbits(2048)
        b = self.rsa.dechunk(self.rsa.chunk(a))
        self.assertEqual(a, b)

    # misc, should not be here
    def test_br(self):
        self.assertEqual(br(1, 8), 128)
        self.assertEqual(br(2, 8), 64)


if __name__ == '__main__':
    unittest.main()

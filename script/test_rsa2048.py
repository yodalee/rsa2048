import unittest
import random
from rsa2048 import *

class TestIntmul2048(unittest.TestCase):
    def test_intmul_2048_result_equals_1(self):
        # Generate random 2048-bit integers a, b, and n
        a = random.getrandbits(2048)
        b = random.getrandbits(2048)
        n = random.getrandbits(2048)

        # Ensure n is odd (for RSA encryption)
        if n % 2 == 0:
            n += 1

        # Calculate the result using the function
        result = intmul_2048(a, b, n)

        # Check if the result is equal to 1
        self.assertEqual(result, (a * b) % n, "Result of intmul_2048 is incorrect")

if __name__ == '__main__':
    unittest.main()
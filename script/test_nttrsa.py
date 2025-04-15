import unittest
from nttrsa import NttRsa


class TestNttRsa(unittest.TestCase):
    def setUp(self):
        # Example setup with N=2048, l=8, len_poly=4, q1=17, q2=19
        self.ntt_rsa = NttRsa(N=2048, l=8, len_poly=4, q1=17, q2=19)

    def test_chunk(self):
        # Test chunking a number
        number = 0x12345678  # Example number
        expected_chunks = [0x78, 0x56, 0x34, 0x12]  # 8-bit chunks
        self.assertEqual(self.ntt_rsa.chunk(number), expected_chunks)

    def test_dechunk(self):
        # Test dechunking a list of integers
        chunks = [0x78, 0x56, 0x34, 0x12]  # Example chunks
        expected_number = 0x12345678  # Original number
        self.assertEqual(self.ntt_rsa.dechunk(chunks), expected_number)

    def test_dechunk_with_overflow(self):
        # Test dechunking with overflow (elements larger than 8 bits)
        chunks = [0x1FF, 0x1FF, 0x1FF, 0x12]  # First 3 elements exceed 8 bits
        expected_number = (0x1FF) + (0x1FF << 8) + (0x1FF << 16) + (0x12 << 24)
        self.assertEqual(self.ntt_rsa.dechunk(chunks), expected_number)

    def test_chunk_dechunk_consistency(self):
        # Test that chunking and dechunking are consistent
        number = 0x12345678  # Example number
        chunks = self.ntt_rsa.chunk(number)
        reconstructed_number = self.ntt_rsa.dechunk(chunks)
        self.assertEqual(reconstructed_number, number)


if __name__ == "__main__":
    unittest.main()

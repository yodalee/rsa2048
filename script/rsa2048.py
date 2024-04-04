
def intmul_2048(a, b, n):
    """Calculate a * b % n with a, b, n are 2048 bits integer"""
    assert(a < 2**2048)
    assert(b < 2**2048)
    assert(n < 2**2048)
    return (a * b) % n


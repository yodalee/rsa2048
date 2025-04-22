# helper script to check the size of a ring
# Below is the note about ring we used for RSA 2048 and 4096 with 32 bits number
# N = 2048; chunking l = 11; poly length n = 384; NTT = 128; q = 12289 x 65537
#  q = 12289 => possible target: 81 ^ 128 = 1 mod 12289
#  q = 65537 => possible target: 2469 ^ 128 = 1 mod 65537
#               another possible target: 4938 ^ 128 = 1 mod 65537
#               since 4936 ^ 4 = 2 mod 65537, this give us extra advantages that
#               most twiddle factors are power of 2
# N = 4096; chunking l = 11; poly length n = 768; NTT = 256; q = 25601 x 65537
#  q = 25601 => possible target: 233 ^ 256 = 1 mod 25601
#  q = 65537 => possible target: 141 ^ 256 = 1 mod 65537
#               another possible target: 5574 ^ 256 = 1 mod 65537
#               same reason that 5574 ^ 8 = 2 mod 65537

def ring(a, r, dump=False):
    if (a <= 1):
        print(f"invalid a = {a}")
        return
    i = 1
    b = a
    while True:
        if dump:
            print(f"{i} = {b}")
        if b == 1:
            break
        i += 1
        b = b * a % r
    return i


if __name__ == "__main__":
    q = 12289
    n = 128
    for i in range(q):
        if ring(i, q) == n:
            print(i)
            break
    ring(i, q, True)

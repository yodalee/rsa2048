

def ring(a, r):
    i = 1
    b = a
    while True:
        print(f"{i} = {b}")
        if b == 1:
            break
        i += 1
        b = b * a % r

if __name__ == "__main__":
    ring(1229, 12289)

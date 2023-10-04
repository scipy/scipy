from random import getrandbits


def gen_random_key(n):
    """Generates a random key of bits (with 0s or 1s) of length n"""
    return bin(getrandbits(n))[2:].zfill(n)

def xor(m, k):
    """Given strings m and k of characters 0 or 1,
    it returns the string representing the XOR
    between each character in the same position.
    This means that m and k should be of the same length.

    Use this function both for encrypting and decrypting!"""
    a = int(m, base=2)
    b = int(k, base=2)
    return bin(a ^ b)[2:].zfill(len(m))

def gen_key(string):
    return "".join(choice(string.printable) for _ in s)

def encrypt(string, key):
    return "".join(chr(ord(i) ^ ord(j)) for (i, j) in zip(string, key))

def decrypt(ciphertext, key):
    return encrypt(ciphertext, key)


if __name__ == "__main__":
    print(xor(
    ))

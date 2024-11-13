import numpy
from dnacryptography import Encryption, Decryption
from simulation import OFDMSystem
import matplotlib.pyplot as plt

text = "Hello, Hi How are you? Good morning"
print(f"[+] Original text: {text}")
encrypted_text, key = Encryption().encrypt(text)
print(f"[+] Encrypted text: {encrypted_text}")
def convert_dna_to_decimal_matrix(string: str, col: int = 8) -> list[list[int]]:
    mapping = {
        "A": 0,
        "T": 1,
        "G": 2,
        "C": 3
    }
    n = len(string)
    row = n // col
    res = []
    c = 0
    for i in range(row):
        temp = []
        for j in range(col):
            temp.append(mapping[string[c]])
            c += 1
        res.append(temp)
    return res

def convert_decimal_matrix_to_dna(matrix: list[list[int]]) -> str:
    reverse_mapping = {
        0: "A",
        1: "T",
        2: "G",
        3: "C"
    }
    dna_string = ""
    for row in matrix:
        for value in row:
            dna_string += reverse_mapping[value]
    return dna_string

matrix = convert_dna_to_decimal_matrix(encrypted_text, col = 4)

print("### Before sending via OFDM convert sequence to decimal ###")

print(matrix)

ofdm = OFDMSystem(4)

matrix = numpy.array(matrix)

transmitted_ofdm_symbols = ofdm.transmit(matrix)

print("### By OFDM transmission function the data is converted to symbols ###")
print(transmitted_ofdm_symbols)

received_bits = ofdm.receive(transmitted_ofdm_symbols)

print("After recieving from the ofdm reciever: ")

print(received_bits)

recieved_cipher = convert_decimal_matrix_to_dna(received_bits)


print(f"Converting decimal matrix to ATCG sequence: {recieved_cipher}")


decrypted_text = Decryption().decrypt(recieved_cipher, key)


print("After decrypttion: ", decrypted_text)


plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
plt.title("Transmitted Data")
plt.imshow(matrix, aspect='auto', cmap='viridis')
plt.subplot(1, 2, 2)
plt.title("Received Data")
plt.imshow(received_bits, aspect='auto', cmap='viridis')
plt.show()
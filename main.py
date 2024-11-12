from dnacryptography import Encryption, Decryption

text = "Hello, Hi How are you?"

print("Original text: ", text)



enc = Encryption()

enc.rounds_no = 10

encrypted_text, key = enc.encrypt(text)

with open("key.txt", "w") as f:
    f.write(key)

print("After encrypttion text: ", encrypted_text)


print("Key: ", key[:200], "....", key[-200:])

print("Length of key: ", len(key))
print("Length of original text: ", len(text))
print("Length of encrypted text: ", len(encrypted_text))

decrypted_text = Decryption().decrypt(encrypted_text, key)

print("After decrypttion: ", decrypted_text)
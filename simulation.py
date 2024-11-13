import numpy as np
import matplotlib.pyplot as plt

class OFDMSystem:
    def __init__(self, modulation_order):
        self.modulation_order = modulation_order

    def transmit(self, data):
        qam_modulated_data = np.exp(1j * 2 * np.pi * data / self.modulation_order)
        ofdm_symbols = np.fft.ifft(qam_modulated_data, axis=1)
        return ofdm_symbols

    def receive(self, ofdm_symbols):
        demodulated_data = np.fft.fft(ofdm_symbols, axis=1)
        received_data = np.round(np.angle(demodulated_data) * self.modulation_order / (2 * np.pi)) % self.modulation_order
        return received_data

num_subcarriers = 64 # row
num_symbols = 100 # column
modulation_order = 4 

#ATGCGETSTSTSTSTSTTYYTYTYTYYTYTYT

"""
A T G C G
G T G T G
T S T S T
S T S T S
T Y Y T Y

"""
# Generate random data
data = np.random.randint(0, modulation_order, (num_symbols, num_subcarriers))

ofdm_system = OFDMSystem(modulation_order)
ofdm_symbols = ofdm_system.transmit(data)
received_symbols = ofdm_symbols
received_data = ofdm_system.receive(received_symbols)

print(data.shape)


# plt.figure(figsize=(12, 6))
# plt.subplot(1, 2, 1)
# plt.title("Transmitted Data")
# plt.imshow(data, aspect='auto', cmap='viridis')
# plt.subplot(1, 2, 2)
# plt.title("Received Data")
# plt.imshow(received_data, aspect='auto', cmap='viridis')
# plt.show()
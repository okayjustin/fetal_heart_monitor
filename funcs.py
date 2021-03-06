import sys
import pyaudio
import wave
import struct
import numpy as np
import RPi.GPIO as GPIO
from scipy import signal
from peakutils.peak import indexes


PLOT = False

if (PLOT):
    import matplotlib.pyplot as plt

#def main(argv):
#    audio_data = input_audio(argv[1])
#    print(audio_data)

# Input: path of CSV file containing h[n] generated by MATLAB
# Output: numpy array of h[n] values
def parse_hn(file_path):
    # hn = np.array([])
    hn_arr = []

    with open(file_path, 'r') as fp:
        file_data = fp.read().strip('\n').split(",")

    for data in file_data:
        hn_arr.append(float(data))

    return hn_arr


# Input: path of .wav audio file
# Output: numpy array of audio data values (only left channel if .wav has more than 1 channel)
def input_audio(file_path):
    audio_data = []

    wave_file = wave.open(file_path, 'r')
    wave_samplewidth = wave_file.getsampwidth()
    wave_chs = wave_file.getnchannels()
    num_frames = wave_file.getnframes()

    for index in range(num_frames):
        wave_data = wave_file.readframes(1)

        # Mono
        if (wave_chs == 1):
            data = struct.unpack("<h", wave_data)
        # Stereo
        elif (wave_chs == 2):
            data = struct.unpack("<hh", wave_data)

        audio_data.append(data[0])

    wave_file.close()

    return np.array(audio_data)


# Input:
#   xn: data sequence to filter
#   hn: impulse response sequence
# Output: filtered sequence
def filter_data(xn, hn):
    yn = np.convolve(xn, hn) #, mode='valid')
    return yn.astype(int)


# Input: xn: data sequence to take short time fourier transform of
# Output: stft result
def do_stft(xn, fs, nperseg):
    noverlap = 0
    nfft = 512  #nperseg;

    s = np.array([])
    f, t, Zxx = signal.stft(xn, fs, nperseg=nperseg, noverlap=noverlap, nfft=nfft)

#    plt.pcolormesh(t, f, np.abs(Zxx), vmin=0) #, vmax=amp)
#    plt.title('STFT Magnitude')
#    plt.ylabel('Frequency [Hz]')
#    plt.xlabel('Time [sec]')
#    plt.show()

    return Zxx

# Input: 2D spectrogram data
# Output: 1D signal power at each time step
def stft2sigpower(s):
    sig_power = np.sum(np.abs(s),0)
    return sig_power

# Input: array of signal powers over time
# Output: peak values in sig_power and locations as indices
def findpeaks(signal, thres=0.2, min_dist=1):
    peaks = np.array([])
    peak_idx = np.array([])

    peak_idx = indexes(signal, thres=thres, min_dist=min_dist)
    peaks = np.array([signal[idx] for idx in peak_idx])
    return peaks, peak_idx

# Encodes a value into an array of 5 LEDs
def encodeVal(val):
    LED_array = [1 if ((val > 105 + 20*i) and (val <= 135 + 20*i)) else 0 for i in range(0,5)]
    if (val < 115):
        LED_array[0] = -1
    elif (val > 205):
        LED_array[4] = -1
    return LED_array

PIN_NUMS = [7,11,13,15,16]

# Initialize RPi GPIO
def initGPIO():
    GPIO.setmode(GPIO.BOARD)
    for pin in PIN_NUMS:
        GPIO.setup(pin, GPIO.OUT)
    for i in range(0, len(PIN_NUMS)):
        GPIO.output(pin, GPIO.HIGH)

# Set LEDs to display a value
def display(val):
    LED_array = encodeVal(val)
    for i in range(0, len(PIN_NUMS)):
        GPIO.output(pin, GPIO.LOW if LED_array[i] else GPIO.HIGH)


class SmoothSeq():
    def __init__(self):
        self.win_len = 5000

        self.summed = 0
        self.buffer = np.zeros(self.win_len)

    # signal: sequence of values to smooth
    def smooth(self, signal):
        smoothed = []

        for i in range(0, self.win_len):
            self.summed = self.summed + signal[i] - self.buffer[i]
            smoothed.append(self.summed)
        for i in range(self.win_len, len(signal)):
            self.summed = self.summed + signal[i] - signal[i - self.win_len]
            smoothed.append(self.summed)

        self.buffer = signal[-self.win_len:]
        return np.array(smoothed) / self.win_len


    def normalize(self, signal):
        norm = np.linalg.norm(signal)
        if norm == 0:
            return signal
        return signal / norm

#if __name__ == "__main__":
#    main(sys.argv)

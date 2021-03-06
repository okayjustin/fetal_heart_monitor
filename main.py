from funcs import *
from recorder import *
import numpy as np
import mute_alsa
import pyaudio
import time
from timeit import default_timer as timer
from helper_funcs import *
import sys, os
import subprocess

SIMULATE = False
SYSTEM_REC = True
PLOT = False

# Fetal heart monitor
class FHM():
    def __init__(self):
        # Import filter h[n] generated by MATLAB
        hn_path = './hn.csv'
        self.hn_arr = parse_hn(hn_path)

        # Recording setup
        self.format = 'S16_LE'
        self.rate = 16000
        # How long to measure before processing
        self.meas_time = 10   # seconds
        # Setup for writing recording to file
        self.wav_name = 'temp'
        self.arecord_args = ['arecord', '', '-f', str(self.format), '-r', str(self.rate), '-d', str(self.meas_time)]

        # Track heart rate
        self.firstMeas = True
        self.beat_time_delta = RunningStat(5)
        self.thresh_low = 100
        self.thresh_high = 220
        self.stdDev_limit = 50
        self.running_win = 5
        self.heart_rate = RunningStat(self.running_win)

    def sysRec(self, buffer_id, block=False):
        wav_path = self.wav_name + str(buffer_id) + '.wav'
        self.arecord_args[1] = wav_path

        if (block):
            os.system('arecord %s -f %s -r %d -d %d -N' % (wav_path, self.format, self.rate, self.meas_time))
        else:
            self.p = subprocess.Popen(self.arecord_args)

    # Determines heart rate then sets LEDs accordingly
    def assessHR(self, audio_data):
        # Determine heart rate
        heart_rate, stdDev = self.processAudio(audio_data)

        # Check rate against thresholds to determine validity
        if ((heart_rate > self.thresh_low) and
            (heart_rate < self.thresh_high) and
            (stdDev < self.stdDev_limit)):

            # If this is the first meas, push multiple times to fill RunningStat
            if (self.firstMeas):
                self.firstMeas = False
                for _ in range(0, self.running_win):
                    self.heart_rate.push(heart_rate)
            else:
                self.heart_rate.push(heart_rate)

            # Display on LEDs
            display(self.heart_rate.winMean())

            print("Heartrate: %0.1f bpm | Inst: %0.1f" % (self.heart_rate.winMean(), heart_rate))
            print("Std. dev: %0.1f bpm\n" % self.heart_rate.winStdDev())
        else:
            print("Out of bounds: %0.1f bpm. Std dev: %0.1f bpm.\n" % (heart_rate, stdDev))

    # FFT based strategy for detecting heart rate
    def processAudio(self, audio_data):
        start = timer()

        # Bandpass filter data
        filtered_data = filter_data(audio_data, self.hn_arr)

        # Smooth data
        nperseg = 64
        s = do_stft(filtered_data, self.rate, nperseg)
        sigpwr = stft2sigpower(s)
#        print(sigpwr)

        # Determine number of peaks in signal power
        peaks, peaks_idx = findpeaks(sigpwr, 0.15, pow(2,12)/nperseg)

        # Measure peak time distance
        self.beat_time_delta.clear()
        for i in range(2, len(peaks_idx)-1):
            self.beat_time_delta.push(peaks_idx[i] - peaks_idx[i-1])
        mean = self.beat_time_delta.mean()
        stdDev = self.beat_time_delta.stdDev()

        # Calculate frequency
        if (mean != 0):
            heart_rate = 60 / (mean * nperseg / self.rate)
        else:
            heart_rate = 0

        end = timer()
#        print("Processing time: %f" % (end-start))

        #Plot everything
        if (PLOT):
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            ax1.plot(sigpwr)
            ax1.scatter(peaks_idx, peaks)
            plt.title('Signal Power vs. Time')
            plt.ylabel('Signal Power')
            plt.xlabel('Time [sec]')
            plt.show()

        return heart_rate, stdDev

    def detectLoop(self):
        buffer_id = 0
        start = timer()
        self.sysRec(buffer_id, block=True)
        end = timer()
        record_time = end - start
        wait_time = 0.05 * record_time

        while (True):
            # Toggle audio buffer
            if (buffer_id == 0):
                buffer_id = 1
                prev_buffer_id = 0
            else:
                buffer_id = 0
                prev_buffer_id = 1

            # Begin recording into next buffer
            self.sysRec(buffer_id)
            time.sleep(wait_time)
            start = timer()

            # Import audio recording
            print("Importing audio... "),
            sys.stdout.flush()
            wav_path = self.wav_name + str(prev_buffer_id) + '.wav'
            audio_data = input_audio(wav_path)
            print("Done.")

            # Processing recording
            self.assessHR(audio_data)

            end = timer()
            print("Processing time: %f" % (end-start))

            while ((end - start) < self.meas_time):
                end = timer()


if __name__ == "__main__":
    print("Fetal Heart Rate Monitor v0.1")
    initGPIO()

    if (PLOT):
        import matplotlib.pyplot as plt

#    for i in range (110, 210):
#        LED_array = encodeVal(i)
#        print("%d: " % i),
#        print(LED_array)

    monitor = FHM()

    if (SIMULATE):
        # Import audio data
        wav_path = './temp0.wav'
        print("Importing audio... "),
        sys.stdout.flush()
        audio_data = input_audio(wav_path)
        print("Done.")

        monitor.assessHR(audio_data)

#        # Pretend to stream data
#        for i in range(0, len(audio_data) / monitor.chunk):
#            audio_data_segment = audio_data[monitor.chunk * i:monitor.chunk * (i+1)-1]
#            monitor.assessHR(audio_data_segment)
    elif (SYSTEM_REC):
        try:
            monitor.detectLoop()
        except KeyboardInterrupt:
            pass


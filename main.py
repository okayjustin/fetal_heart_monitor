from funcs import *
from recorder import *
import numpy as np
import mute_alsa
import pyaudio
import time
from timeit import default_timer as timer
from helper_funcs import *
import sys, os

SIMULATE = False
SYSTEM_REC = True
PLOT = False

# Fetal heart monitor
class FHM():
    def __init__(self):
        # How long to measure before processing chunk
        self.measwindow = 10   # seconds

        # Import filter h[n] generated by MATLAB
        hn_path = './hn.csv'
        self.hn_arr = parse_hn(hn_path)

        # Recording setup
##        self.p = pyaudio.PyAudio()
#        info = self.p.get_host_api_info_by_index(0)
#        numdevices = info.get('deviceCount')
#        for i in range(0, numdevices):
#                if (self.p.get_device_info_by_host_api_device_index(0, i).get('maxInputChannels')) > 0:
#                    print "Input Device id ", i, " - ", self.p.get_device_info_by_host_api_device_index(0, i).get('name')
#        for i in range(0, numdevices):
#                if (self.p.get_device_info_by_host_api_device_index(0, i).get('maxOutputChannels')) > 0:
#                    print "Output Device id ", i, " - ", self.p.get_device_info_by_host_api_device_index(0, i).get('name')
        self.format = pyaudio.paInt16
        self.rate = 16000 #44100
        self.chunk = int(self.rate * self.measwindow)
        self.meas_time = float(self.chunk) / self.rate

        # Setup for writing recording to file
        self.fname = 'test.wav'
        self.mode = 'wb'
        self.channels = 1

        # Init smoother
#        self.smoother = SmoothSeq()

        # Track heart rate
        self.firstMeas = True
        self.beat_time_delta = RunningStat(5)
        self.thresh_low = 100
        self.thresh_high = 220
        self.stdDev_limit = 50
        self.running_win = 5
        self.heart_rate = RunningStat(self.running_win)

    # Start recording call back
    def start_recording(self):

    #    self.wavefile = self._prepare_file(self.fname, self.mode)
        self.stream = self.p.open(format=self.format,
                        channels= self.channels,
                        rate=self.rate,
                        input=True,
                        output=False,
                        input_device_index=2,
                        output_device_index=1,
                        frames_per_buffer=self.chunk,
                        stream_callback=self.get_callback())
        self.last_cb_time = timer()
        print("Starting recording... ")

    def get_callback(self):
        def callback(in_data, frame_count, time_info, status):
            # Measure time between callbacks
            current_time = timer()
            print("Callback delta: %0.1f" % (current_time - self.last_cb_time))
            self.last_cb_time = current_time

            start_time = timer()
#            audio_data = np.fromstring(in_data, dtype=np.int16)
#            self.assessHR(audio_data)

 #           print("Writing to file..."),
 #           sys.stdout.flush()
 #           self.wavefile.writeframes(in_data)
 #           print("Done.")

            end_time = timer()
            print("Callback %0.1f" % ((end_time - start_time)))
            return in_data, pyaudio.paContinue
        return callback

    def close(self):
        try:
            self.stream.close()
        except:
            pass
        self.p.terminate()

    def sysRec(self):
        wav_path = './temp.wav'
        os.system('arecord %s -f S16_LE -r %d -d %d' % (wav_path, self.rate, self.measwindow))
        return wav_path

    def _prepare_file(self, fname, mode='wb'):
        wavefile = wave.open(fname, mode)
        wavefile.setnchannels(self.channels)
        wavefile.setsampwidth(self.p.get_sample_size(pyaudio.paInt16))
        wavefile.setframerate(self.rate)
        return wavefile

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

            print("Heartrate: %0.1f bpm | Inst: %0.1f" % (self.heart_rate.winMean(), heart_rate))
            print("Std. dev: %0.1f bpm\n" % self.heart_rate.winStdDev())
        else:
            print("Out of bounds: %0.1f bpm. Std dev: %0.1f bpm.\n" % (heart_rate, stdDev))

    # Moving average strategy for detecting heart rate
    def processAudioMA(self, audio_data):
        start = timer()

        # Calculate length of audio file
        total_time = 1.0 * np.size(audio_data) / self.rate

        # Bandpass filter data
        filtered_data = filter_data(audio_data, self.hn_arr)

        # Get data magnitude
        abs_filt_audio_data = abs(filtered_data)
        abs_audio_data = abs(audio_data)

        # Smooth data
        smoothed_abs_filt_audio_data = self.smoother.smooth(abs_filt_audio_data)
        smoothed_abs_audio_data = self.smoother.smooth(abs_audio_data)

        # Determine number of peaks in signal power
        peaks_filt, peaks_idx_filt = findpeaks(smoothed_abs_filt_audio_data, 0.2, 10000)
        peaks, peaks_idx = findpeaks(smoothed_abs_audio_data, 0.2, 10000)

        # Calculate frequency with filter
        num_peaks = np.size(peaks_filt)
        freq = num_peaks / total_time
        print("Heart rate: %0.1f bpm" % (freq*60))

        # Calc freq w/o filter
        num_peaks = np.size(peaks)
        freq = num_peaks / total_time
        print("Heart rate (no filter): %0.1f bpm" % (freq*60))

        end = timer()
        print(end-start)

        #Plot everything
        fig = plt.figure()
        ax1 = fig.add_subplot(311)
        ax1.plot(abs_filt_audio_data, 'b')
        ax1.plot(smoothed_abs_filt_audio_data, 'r')
        ax1.plot(smoothed_abs_audio_data, 'g')

        ax2 = fig.add_subplot(312)
        ax2.plot(smoothed_abs_filt_audio_data, 'r')
        ax2.scatter(peaks_idx_filt, peaks_filt)

        ax3 = fig.add_subplot(313)
        ax3.plot(smoothed_abs_audio_data, 'g')
        ax3.scatter(peaks_idx, peaks)

        plt.title('Signal Power vs. Time')
        plt.ylabel('Signal Power')
        plt.xlabel('Time [sec]')
        plt.show()


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
        peaks, peaks_idx = findpeaks(sigpwr, 0.15, pow(2,13)/nperseg)

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
#        print("Heart rate: %0.1f bpm | Std dev: %0.1f bpm.\n" % (heart_rate, stdDev))

        end = timer()
#        print("Processing time: %f" % (end-start))

        #Plot everything
#        fig = plt.figure()
#        ax1 = fig.add_subplot(111)
#        ax1.plot(sigpwr)
#        ax1.scatter(peaks_idx, peaks)
#        plt.title('Signal Power vs. Time')
#        plt.ylabel('Signal Power')
#        plt.xlabel('Time [sec]')
#        plt.show()

        return heart_rate, stdDev

if __name__ == "__main__":
    print("Fetal Heart Rate Monitor v0.1")

    if (PLOT):
        import matplotlib.pyplot as plt

#    for i in range (110, 210):
#        LED_array = encodeVal(i)
#        print("%d: " % i),
#        print(LED_array)

    monitor = FHM()

    if (SIMULATE):
        # Import audio data
        wav_path = './temp.wav'
        print("Importing audio... "),
        sys.stdout.flush()
        audio_data = input_audio(wav_path)
        print("Done.")

        # Pretend to stream data
        for i in range(0, len(audio_data) / monitor.chunk):
            audio_data_segment = audio_data[monitor.chunk * i:monitor.chunk * (i+1)-1]
            monitor.assessHR(audio_data_segment)
    elif (SYSTEM_REC):
        wav_path = monitor.sysRec()
        print("Importing audio... "),
        sys.stdout.flush()
        audio_data = input_audio(wav_path)
        print("Done.")
        monitor.assessHR(audio_data)
    else:
        # Record until quit
        try:
            monitor.start_recording()
            while(True):
                pass
        except KeyboardInterrupt:
            pass

#    monitor.close()


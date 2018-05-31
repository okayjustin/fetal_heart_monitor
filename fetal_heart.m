clear all
close all
clc
figure_num = 1;

% Import audio file
[Y, fsample] = audioread('fetal_heart.wav');

% Optimal FIR design ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Config parameters
As = 60; % dB
Ap = 0.5; % dB

% Calculate passband ripple and stopband attenuation
delta_s = 10^(-As/20);
delta_p = (10^(Ap/20)-1)/(10^(Ap/20)+1);

% Parameters for optimal filter design
F = [3500 4000 6000 6500];
A = [0 1 0];
DEV = [delta_s delta_p delta_s];
[N,Fo,Ao,W] = firpmord(F,A,DEV,fsample)
N = N + 17; % Increase order a bit

% Optimal filter design
Bk = firpm(N,Fo,Ao,W);

% Analyze filter
M = N + 1
Ak = zeros(1,M);
Ak(1) = 1;
num_of_f_points = fsample / 2;
num_of_n_points = M;
[poles,zeros,HF,Fd,hn,n,figure_num] = show_filter_responses(Ak,Bk,fsample,num_of_f_points,num_of_n_points,figure_num);
csvwrite('hn.csv',hn)
% Filter the audio to remove low freq noise ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Y_left = Y(1:end,1);
Y_right = Y(1:end,2);
% Y_left = Y_left(500000:941000);

% Bandpass filter the signal
Y_left_filtered = filter(Bk, Ak, Y_left);

% Spectrogram settings
win_size = 4000;
fft_len = 64;

% Plot the spectrogram before filter
figure(figure_num)
figure_num = figure_num + 1;
colormap gray
spectrogram(Y_left,win_size,0,fft_len,fsample);
S_nofilt = spectrogram(Y_left,win_size,0,fft_len,fsample);
% win_len = win_size / fsample;
% t = 0:win_len:size(Y_left_filtered,1)/fsample - win_len;
% t_total = size(Y_left_filtered,1)/fsample;

% Plot the spectrogram after filter
figure(figure_num)
figure_num = figure_num + 1;
colormap gray
spectrogram(Y_left_filtered,win_size,0,fft_len,fsample);
S = spectrogram(Y_left_filtered,win_size,0,fft_len,fsample);
win_len = win_size / fsample;
t = 0:win_len:size(Y_left_filtered,1)/fsample - win_len;
t_total = size(Y_left_filtered,1)/fsample;

% Save filtered audio
%audiowrite('fetal_heart_filtered.wav', Y_left_filtered, fsample);

% Sum the magnitudes along the frequency axis to produce a single power at
% each time
summed_nofilt = sum(abs(S_nofilt),1) / fft_len;
summed = sum(abs(S),1) / fft_len;

% Unfiltered ~~~~~~~~~~~
% Find peaks in the waveform
[peaks_nofilt, location_nofilt] = findpeaks(summed_nofilt, 'MinPeakProminence', 0.01);

% Plot the signal power over time
figure(figure_num)
figure_num = figure_num + 1;
plot(t,summed_nofilt)
xlabel('Time (s)');
ylabel('Signal Power dB/Hz');
title('Signal Power vs Time (Unfiltered');

% Overlay the scatterplot of peaks
hold on;
scatter((location_nofilt-1) * win_len, peaks_nofilt);

% Calculate frequency
num_peaks_nofilt = size(peaks_nofilt,2);
heart_freq_nofilt = num_peaks_nofilt / t_total;
heart_bpm_nofilt = heart_freq_nofilt * 60

% Filtered~~~~~~~~~~~
% Find peaks in the waveform
[peaks, location] = findpeaks(summed, 'MinPeakProminence', 0.01);

% Plot the signal power over time
figure(figure_num)
figure_num = figure_num + 1;
plot(t,summed)
xlabel('Time (s)');
ylabel('Signal Power dB/Hz');
title('Signal Power vs Time');

% Overlay the scatterplot of peaks
hold on;
scatter((location-1) * win_len, peaks);

% Calculate frequency
num_peaks = size(peaks,2);
heart_freq = num_peaks / t_total;
heart_bpm = heart_freq * 60

close all; clear all;

data = readmatrix('S1_21_Male.csv', 'Delimiter', '\t');
data = data(:,4);
t = 1:length(data);
acquired_data = [];
thenvelopewindow=128*20;
x = highpass(data, 0.5, 125);
x = lowpass(x, 35, 125);
[rms_signal,maxChunkTime] = envelope(x, thenvelopewindow, 'rms');

%th = triangleThreshold(rms_signal, 24);
%action = (rms_signal>th).*max(rms_signal);
t_1 = 1:length(rms_signal);
plot(t,data)
figure
plot(t,x)
figure
plot(t_1, rms_signal)

Fs = 125; % Sampling frequency
T = 1/Fs; % Sampling period
L = length(x); % Length of signal
t = (0:L-1)*T; % Time vector
Y = fft(x); % Fourier transform
P2 = abs(Y/L); % Two-sided spectrum
P1 = P2(1:L/2+1); % Single-sided spectrum
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L; % Frequency vector
figure
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

%Question 9
%Neelabhro Roy
clear all;
close all;
clc;

clear y Fs;
%Here we are reading the Audio signal
%[y,Fs] = audioread('DITTY1.WAV');
[y,Fs] = audioread('BUMMER.WAV');
% This returns the sampled data into y, and the sampling rate of the data
% to Fs

sound(y,Fs);
hold on

Y = ['The Sample rate of the signal is ',num2str(Fs)];   
disp(Y);
%plot(length(y)./1000,y);

specgram(y,length(y),Fs);
title('Spectrogram of the Input Audio Signal');
figure;
%The  Spectrogram of the Audio signal

plot(psd(spectrum.periodogram,y,'Fs',Fs,'NFFT',length(y)));
% This graph provides us with the Power Spectral Density of the audio
% signal

figure;

N = length(y);
FFT = fft(y) ./ N;
%Normalising the Fourier Transform
Fn = Fs/2; 
% The Nyquist Frequency
Freq = (( linspace(0,1,fix(N/2)+1)) .* Fn);
Index = 1 : length(Freq);
%Assigning the indices
stem(Freq./1000, abs(FFT(Index))*2);

title('Fourier Transform of the Audio signal');
xlabel('Frequency (kHz)');
ylabel('Amplitude');
figure;

plot(20*log10(abs(FFT)));
title('Graph for Bandwidth of the Signal');
%xlabel('Frequency (kHz)');
%ylabel('Amplitude');
figure;
%For bandwidth of the signals
obw(y,Fs);

%sound(real(y),Fs);
figure;
%d1 = fdesign.lowpass('Fp,Fst,Ap,Ast',0.0001,0.1,1,60);


d1 = fdesign.lowpass('N,Fc',10,1200,48000);
% We are performing Lowpass filtering at the Half Bandwidth
designmethods(d1);

f1 = design(d1, 'window');  
%fvtool(f);
Q = filter(f1,y);

stem(Freq./1000,abs(Q(Index)));
figure;

obw(Q,Fs);
figure;


load handel.mat

filename = 'handel1.wav';
audiowrite(filename,Q,Fs);
clear Q Fs

[Q,Fs] = audioread(filename);
sound(Q,Fs);

hold on

specgram(Q,length(Q),Fs);
title('Spectrogram of the Lowpass Filtered Audio Signal');
figure;



d2 = fdesign.highpass('N,Fc',10,1200,48000);
%('Fp,Fst,Ap,Ast',0.0001,0.1,1,60);
% We are doing Highpass Window filtering of the input signal, and as is
% evident from the resultant Sound, only the High frequencies above the
% Half Bandwidth are resolved and audible.


designmethods(d2);
f2 = design(d2, 'window');  
%fvtool(f);

I = filter(f2,y);
stem(Freq./1000,abs(I(Index)));
figure;

obw(I,Fs);
figure;

load handel.mat

filename = 'handel2.wav';
audiowrite(filename,I,Fs);
clear I Fs

[I,Fs] = audioread(filename);
sound(I,Fs);

specgram(I,length(I),Fs);
title('Spectrogram of the Highpass Filtered Audio Signal');
clc
clear
close all
 
%play a single tone
Fs = 8000;         % Sampling rate Fs:  This must be at least 2x your highest frequency!!!
Ts = 1/Fs;             % Compute sampling time Ts
N  = 8000;            % 8000 number of time samples separated by Ts seconds

%mary had a little lamb
t = Ts*[0:1:N/2-1];            % Form t vector with values separated by Ts seconds
t2 = Ts*[0:1:N/2.2-1];
t3 = Ts*[0:1:1.7*N-1];
tP = Ts*[0:1:N/10-1];
f0 = 1000;                   % fundamental frequency
soundE = 1*sin(2*pi*f0*t); % Compute the samples of the sinewave
soundE2 = 1*sin(2*pi*f0*t2); 
soundD = 1*sin(2*pi*(0.9*f0)*t);
soundD2 = 1*sin(2*pi*(0.9*f0)*t2);
soundC = 1*sin(2*pi*(0.8*f0)*t); 
soundC2 = 1*sin(2*pi*(0.8*f0)*t3); 
pause = 1*sin(2*pi*0*f0*tP); 
bigpause = 1*sin(2*pi*0*f0*t);
soundA = 1*sin(2*pi*(1.2*f0)*t2);  
song = [soundE,soundD,soundC,soundD,soundE2,pause,soundE2,pause,soundE2,bigpause,soundD2,pause,soundD2,pause,soundD2,bigpause,soundE2,pause,soundA,pause,soundA,bigpause,soundE,soundD,soundC,soundD,pause,soundE2,pause,soundE2,pause,soundE2,pause,soundE,pause,soundD,pause,soundD,pause,soundE,pause,soundD,pause,soundC2];
 
sound(song, Fs) % Play the song 

%%%% write song to a wav file  %%%
audiowrite('Mary Had A Little Lamb.wav',song,Fs)

%don't play chord until song is over
y=audioplayer(song,Fs);
playblocking(y);  
 
%play a single tone
Fs = 8000;         % Sampling rate Fs:  This must be at least 2x your highest frequency!!!
Ts = 1/Fs;             % Compute sampling time Ts
N  = 8000;            % 8000 number of time samples separated by Ts seconds
 
t = Ts*[0:1:N-1];            % Form t vector with values separated by Ts seconds
f0 = 1000;                   % fundamental frequency
soundt = 1*sin(2*pi*f0*t);       % Compute the samples of the sinewave
 
%%%  read in a wav file (Dchord) %%
[soundt2, Fs] = audioread('chord.wav');

%compute the  frequency domain signal
NFFT = 2^nextpow2(length(soundt2)); % Next power of 2 from length of y - this will set your number of FFT points. 
spectrum = fft(soundt2,NFFT)/length(soundt2);  %take the FFT to compute the spectrum
f = Fs/2*linspace(0,1,NFFT/2+1);     %determine the frequencies within the spectrum.  This will compute the correct frequencies
 
%plot the spectrum of the tone over frequency
figure
plot(f/1e3,2*abs(spectrum(1:NFFT/2+1)))
title('Single-Sided Amplitude Spectrum of Chord')
xlabel('Frequency (kHz)')
ylabel('|Y(f)|')
axis([0 5 0 .03])
 
Ts=1/Fs;
Nt=length(soundt2);
chordt=Ts*[0:1:Nt-1];
 
%plot the spectrum of the tone over time
figure
plot(chordt,soundt2)
title('Single-Sided Amplitude Spectrum of Chord')
xlabel('Time (seconds)')
ylabel('|Y(f)|')
axis([0 10 0 .8])
 
%play the chord
sound(soundt2, Fs)
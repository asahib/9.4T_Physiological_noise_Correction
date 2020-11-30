clear all
close all

%load physiological data
load puls
load respr
pulse= puls;
resp=respr;

%Physiological data - Fourier analysis-------------------------------------
Fs=50;
T=1/Fs;

% t_resp=(0:L_resp-1)*T;
% t_pulse=(0:L_pulse-1)*T;

nfft_pulse=2^nextpow2(length(pulse));
nfft_resp=2^nextpow2(length(resp));

f_pulse= 0:Fs/length(pulse):25;
f_resp= 0:Fs/length(resp):25;

Fpulse=fft(pulse,nfft_pulse)/length(pulse);
Fresp=fft(resp,nfft_resp)/length(resp);

abs_Fpulse=abs(Fpulse(1:nfft_pulse/2+1));
abs_Fresp=abs(Fresp(1:nfft_resp/2+1));

plot(f_pulse,2*abs_Fpulse(1:nfft_pulse/2+1))
title('Single-Sided Amplitude Spectrum of respiration')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')

% figure(1)
% subplot(1,2,1)
% plot(f_resp,2*abs_Fpulse(1:nfft_pulse/2+1))
% xlabel('Frequency (Hz)')
% ylabel('|Y(f)| (Pulse)')
% title('Single-Sided Amplitude Spectrum of the Pulse data')
% subplot(1,2,2)
% plot(f,2*abs_Fresp(1:nfft_resp/2+1))
% xlabel('Frequency (Hz)')
% ylabel('|Y(f)| (Respiration)')
% %text(-10,10.2,'Test title spanning two subplots -- Some fine tuning will be required')
% title('Single-Sided Amplitude Spectrum of the Respiration data')


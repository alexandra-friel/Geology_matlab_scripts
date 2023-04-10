%%
%Name: Alexandra Friel 
%Date: March 6th, 2023
%Class: Digital Signal Processing
%This script aims to Understand how to use and implement MATLAB’s fft
%function and introduce the concept of a spectrogram.

%%
% Given vector: 
x = [0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0];

X = dft(x);
phase_spect = angle(X);
amp_spect = abs(X);

freq = linspace(0, pi, (length(x)/2) + 1);
K = length(x)/2 + 1;


subplot(2,1,1)
plot(freq, amp_spect,'k','LineWidth',2)
xlabel('Frequency (rad)')
ylabel('Amplitude')
title('Amplitude Spectrum')
grid on; 

subplot(2,1,2)
plot(freq, phase_spect,'k','LineWidth',2)
xlabel('Frequency (radians)')
ylabel('Phase')
title('Phase Spectrum')
grid on; 

%%
x = [0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0];
X = fft(x);

N = length(x);
freq_x = linspace(0, pi, floor(N/2)+1);
amp_spect = abs(X(1:N/2+1));


subplot(2,1,1)
plot(freq_x, amp_spect, 'k','LineWidth',2);
xlabel('Frequency (rad)');
ylabel('Amplitude');
title('Amplitude Spectrum - FFT');
grid on; 

subplot(2,1,2)
phase_spect = angle(X(1:N/2+1));
plot(freq_x, phase_spect,'k','LineWidth',2)
xlabel('Frequency (rad)');
ylabel('Phase');
title('Phase Spectrum -  FFT');
grid on; 

%%

%dft will calculate and return a complex vector that includes both the magnitude and phase of the Fourier coefficients.
%fft will compute a complex vector that only includes the magnitude of the Fourier coefficients,
%and not  the phase. 


%% 
% Define data vector x
x = [0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0];

n = length(x);
x_padded = [x zeros(1, mod(n, 2))];
tic;

X_dft = fft(x_padded);
time_dft = toc;
tic;

X_fft = fft(x_padded);
time_fft = toc;

fprintf('time in s for DFT: %f seconds', time_dft);
fprintf('time in s for FFT: %f seconds', time_fft);

%% 
load('whistle.mat')
y = Y; 

N = length(y);
Fs = 44100;
steps= Fs/N;

t = (0:N-1) * steps;
Y = fft(y); 
Y_amps = abs(Y); 

figure;
subplot(2,1,1)
plot(t, Y, 'k')
xlim([0 10000]) %Cropped freq domain
ylabel('Amplitude')
title('Time-domain signal')

subplot(2,1,2)
plot(t, Y_amps, 'k')
xlim([0 10000])
xlabel('Time (s)')
ylabel('Magnitude')
title('Frequency-domain spectrum')

%%
spectrogram(Y_amps(1:length(Y_amps)/2),256,50,1024)
title('Spectrogram from whistle data')
colormap(jet); %I was able to see more variations with a different color map.

%This did not appear to be stationary.


%This was so cool! Thank you!

%% Function below, do not alter!


function X = dft(x)
% compute the DFT of a data vector x using basis functions
if size(x,1)<size(x,2)
x = x'; % transpose matrix if needed
end
K = length(x)/2 + 1; % number of frequency points to compute
fax = linspace(0,pi,K); % frequency axis in radians
n = (0:length(x)-1)'; % vector of time series indices (starting with 0)

for k=1:K % loop through each frequency (in radians)
ff = fax(k); % incremental frequency
c(:,k) = cos(ff*n); % cosine basis function
s(:,k) = sin(ff*n); % sine basis function
re(k) = sum(x.*c(:,k)); % real coefficient
im(k) = sum(x.*s(:,k)); % imaginary coefficient
end
X = re - im*(i); % complex spectrum
end 
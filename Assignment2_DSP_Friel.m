%% 
%Alexandra Friel 
%DSP assignment 2 


%% 
freq = 1.9;
N = 10000;
sps = 1000;
time= linspace(0, (N-1)/sps, N);
w = freq*2*pi;
x = time*w;

subplot(4,1,1) 
y = sin(x);
plot(time,y,'m.','LineWidth',1)
xlabel('time (s)')
Nyquist1 = sps/2;
title(['New nyquist frequency = ' num2str(Nyquist1)])


subplot(4,1,2)
dec_factor = 50;
x = linspace(0, 19*2*pi, 10000);
y = sin(x);
t_down = time(1:dec_factor:end);
y_down = y(1:dec_factor:end);

plot(t_down, y_down, '-m')
xlabel('Time (s)')
nyquist_freq_down = sps/(2*dec_factor);
title(['New Nyquist Frequency = ' num2str(nyquist_freq_down)])

subplot(4,1,3)
dec_factor = 500;
x = linspace(0, 19*2*pi, 10000);
y = sin(x);
t_down = time(1:dec_factor:end);
y_down = y(1:dec_factor:end);

plot(t_down, y_down, '-m')
xlabel('Time (s)')
nyquist_freq_down = sps/(2*dec_factor);
title(['New Nyquist Frequency = ' num2str(nyquist_freq_down)])


subplot(4,1,4)
dec_factor = 476;
x = linspace(0, 19*2*pi, 10000);
y = sin(x);
t_down = time(1:dec_factor:end);
y_down = y(1:dec_factor:end);

plot(t_down, y_down, '-m')
xlabel('Time (s)')
nyquist_freq_down = sps/(2*dec_factor);
title(['New Nyquist Frequency = ' num2str(nyquist_freq_down)])

%%
load whistle

%% 
ts = (1:length(Y))/Fs;
subplot(3,2,1)
plot(ts,Y)
%sound(Y,Fs)

title('Original Whistle')
xlabel('time (s)')
Fs_nyquist = Fs/2;

subplot(3,2,2) 
n = 2^nextpow2(length(Y));
Y_fft = fft(Y,n);
f = (0:n-1)*(Fs/2)/n;
plot(f/1000,abs(Y_fft))
xlabel('frequency (kHz)')
hold on

Nyquist = Fs/2;
plot([Nyquist Nyquist]/1000, ylim, 'm-')
plot([Nyquist Nyquist]/(1000*12), ylim, 'm-')
title('Frequency Spectrum (Original)')

subplot(3,2,3)
res = 12;
Yres = downsample(Y, res);
tres = downsample(ts, res);

%sound(Yres,Fs/res);
plot(tres,Yres, 'Color','b');
title('downsampled whistle time series')
xlabel('time (s)')

subplot(3,2,4) 
nres = 2^nextpow2(length(Yres));
fsres = (0:nres-1)*(Fs/res)/nres;
Ares = abs(fft(Yres, nres));
plot(fsres/1000, Ares,'Color',[128, 0, 0]/255);
xlabel('Frequency (kHz)')
ylabel('Amplitude')

xlabel('frequency (kHz)')
hold on
Fnres = Fs/res/2;

plot([Fnres Fnres]/1000,ylim,'m');
title('downsampled frequency spectrum')
axis tight

subplot(3,2,5)
dec = 12;
Ydec = decimate(Y,dec);
tdec = ts(1:dec:end);
%sound(Ydec,Fs/dec);

plot(tdec,Ydec)
axis tight
title('decimated whistle time series')
xlabel('time (s)')

subplot(3,2,6) 
ndec = 2^nextpow2(length(Ydec));
Adec = abs(fft(Ydec,ndec));
fsdec = ((0:ndec-1)/ndec)*Fs/dec;
plot(fsdec/1000,(Adec),'Color',[128, 0, 0]/255);

xlabel('frequency (kHz)')
hold on
Fndec = Fs/dec/2;

plot([Fndec Fndec]/1000,ylim,'m-')
title('decimated frequency spectrum')


%%
subplot(1, 1, 1)
plot(fsres/1000, Ares*12, 'Color',[128, 0, 0]/255);
hold on
plot(fs/1000, A, 'b', 'LineWidth', 2, 'Color', [1, 0.5, 0]);
xlabel('Frequency (kHz)')
ylabel('Scaled Amplitude')
xlim([0, 1.5*Fnres/1000])

title('Original and Aliased Frequency Spectra')
legend('Aliased Whistle', 'Original Whistle', 'Location', 'northwest')
line([Fnres/1000, Fnres/1000], ylim, 'Color', 'k', 'LineWidth', 1)


%%

subplot(2, 2, 1)
plot(ts, Y,'Color',[128, 0, 0]/255);
xlabel('Time (s)')
title('Original Signal')
subplot(2, 2, 2)
[B, F, T] = spectrogram(Y, 1024, 512, 1024, Fs);
imagesc(T, F/1000, log10(abs(B)))
set(gca, 'YDir', 'normal')
ylabel('Frequency (kHz)')
title('Spectrogram Original Signal')

subplot(2, 2, 3)
plot(tres, Yres,'Color',[128, 0, 0]/255);
xlabel('Time (s)')
title('Downsampled Signal')

subplot(2, 2, 4)
[B, F, T] = spectrogram(Yres, 256, 128, 1024, Fs/res);
imagesc(T, F/1000, log10(abs(B)))
set(gca, 'YDir', 'normal')
ylabel('Frequency (kHz)')
title('Spectrogram of Downsampled Signal')




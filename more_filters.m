load('Erebus_seismogram.mat')

%%

figure(1); clf
sps = hdr.sps
ts = (1:length(data))/sps;
sens = 3200; % V/m/s
atod = hdr.atod; % V/count
vel = detrend(data,'constant')*hdr.atod/sens*1e6; % convert to micrometers/s
plot(ts,vel,'k')
ylabel('\mum/s')
xlabel('time (s)')
title('Erebus velocity seismogram')

%% filter seismogram using Butterworth type filter

npoles = 2;
fc = .25; % Hz low-pass filter - what period is this?
dfc = fc/(sps/2);
[B,A] = butter(npoles,dfc,'low')
vfilt1 = filtfilt(B,A,vel);
hold on
plot(ts,vfilt1,'linewidth',3)
axis tight

%% filter seismogram using running average filter with 4 second avg

h = ones(4*sps,1)/(4*sps);
vfilt2 = conv(vel,h,'same')
hold on
plot(ts,vfilt2,'linewidth',3)
axis tight
legend('raw','lowpass','running average')

%% plot frequency response of Butterworth filter using coefficients

figure(2); clf
[H,f] = freqz(B,A,1024)
fHz = f/pi*(sps/2); % convert to Hz
subplot(2,1,1)
semilogx(fHz,abs(H))
grid on
xlim([.01 20])
ylabel('amplitude response')
title('Bode plot for filter')
xlabel('frequency (Hz)')

subplot(2,1,2)
semilogx(fHz,angle(H)*180/pi)
grid on
xlim([.01 20])
xlabel('frequency (Hz)')
ylabel('pahse (degrees)')

%% next step is to plot frequency spectrum for original data and filtered data

figure(3); clf
VEL = fft(vel);
semilogx(abs(VEL))
hold on
VFILT1 = fft(vfilt1);
semilogx(abs(VFILT1))
VFILT2 = fft(vfilt2);
semilogx(abs(VFILT2))
grid on
legend('raw','lowpass','running average')

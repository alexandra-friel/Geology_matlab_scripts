

%% Load Wind Speed Data & design filters

% 1 Hr = 3600 s
% ts = 1: length(data) 
%plot(ts,data)
%ts_hr = ts /24
%plot(ts_hr, data) make sure axis is properly scaled
%sample rate = 1/delta t -> 1/hr = *24hrs/day = 24/day

sps = 24; % samples per second
time = (1:length(HzWindSpd_T))/sps;
time_24hrs = time/24; 
wind_speed_detrended = detrend(HzWindSpd_T, 'constant'); % remove mean

low_pass_cutoff = 0.25; % Hz
high_pass_cutoff = 0.25; % Hz
num_poles = 2;

[lowpassB, lowpassA] = butter(num_poles, low_pass_cutoff/(sps/2), 'low'); % low-pass Butterworth filter
[highpassB, highpassA] = butter(num_poles, high_pass_cutoff/(sps/2), 'high'); % high-pass Butterworth filter
[bandpassB, bandpassA] = butter(num_poles, [0.1/(sps/2), 5/(sps/2)], 'bandpass'); % band-pass Butterworth filter

wind_speed_lowpass = filtfilt(lowpassB, lowpassA, wind_speed_detrended);
wind_speed_highpass = filtfilt(highpassB, highpassA, wind_speed_detrended);
wind_speed_bandpass = filtfilt(bandpassB, bandpassA, wind_speed_detrended);

%% Plot Filtered Data

figure(2); 
% Low-pass filter

plot(time, wind_speed_detrended, 'k')
hold on
plot(time, wind_speed_lowpass, 'linewidth', 2)
ylabel('Wind Speed (mph)')
xlabel('Time in days')
title('Low-pass Filtered Wind Speed')
axis tight
xlim([0 580])
legend('Original Data', 'Low-pass Filtered')

%%
figure(2); 
% Low-pass filter

plot(time, wind_speed_detrended, 'k')
hold on
plot(time, wind_speed_lowpass, 'linewidth', 2)
ylabel('Wind Speed (mph)')
xlabel('Time in days')
title('Low-pass Filtered Wind Speed')
axis tight

% Stretch out the x-axis by increasing the range
xlim([time(1)-10, time(end)+10])

legend('Original Data', 'Low-pass Filtered')



% subplot(2,1,2)
% plot(time_24hrs,wind_speed_detrended,'k','LineWidth',1)
% ylabel('Wind Speed (m/s)')
% xlabel('Time in days')
% title('Wind Speed from Lower Deer Point, Idaho')
% legend('Original Data')


% % High-pass filter
% subplot(3,1,2)
% plot(time, wind_speed_detrended, 'k')
% hold on
% plot(time, wind_speed_highpass, 'linewidth', 1)
% ylabel('Wind Speed (mph)')
% xlabel('Time (s)')
% title('High-pass Filtered Wind Speed')
% axis tight
% legend('Original Data', 'High-pass Filtered')
% 
% % Band-pass filter
% subplot(3,1,3)
% plot(time, wind_speed_detrended, 'k')
% hold on
% plot(time, wind_speed_bandpass, 'linewidth', 1)
% ylabel('Wind Speed (mph)')
% xlabel('Time (s)')
% title('Band-pass Filtered Wind Speed')
% axis tight
% legend('Original Data', 'Band-pass Filtered')

%% Plot Frequency Response of Low-pass Filter - from  Jeff

t0 = 6.655; % time scale to decay to
fc = 1/2/pi/t0; % definition of digital corner
npoles = 4;
sps = 24; % suppose this is the sample rate then...
corner = fc*sps;

Wn = fc*2 % Wn is fraction of Nyquist whereas Fc is digital frequency
[Blow,Alow] = butter(npoles,Wn,'low')
[Bhigh,Ahigh] = butter(npoles,Wn,'high')
figure(1); clf

[H,f] = freqz(Blow,Alow,1024)
fHz = f/pi/2*sps
subplot(3,1,1) % amplitude response
semilogx(fHz,abs(H))

grid on
xlim([.01 sps/2])
ylabel('amplitude response')
title('Bode plot for filter')
xlabel('frequency (Hz)')
hold on

plot(xlim,sqrt(2)/2+[0 0])
plot(corner+[0 0],ylim)
ylim([0 1])
subplot(3,1,2) % phase response
semilogx(fHz,angle(H)*180/pi)
grid on

xlim([.01 sps/2])
xlabel('frequency (Hz)')
ylabel('phase (degrees)')
subplot(3,1,3) % convert amplitude response to dB
dB = 20*log10(abs(H))
semilogx(fHz,dB)
grid on

xlim([.01 sps/2])
ylabel('dB')
xlabel('frequency (Hz)')
hold on
plot(xlim,[-3 -3])
plot(corner+[0 0],ylim)
ylim([-30 0])
set(gca,'ytick',[-30:6:0])

%% Impulse Response for wind data: 

low_pass_cutoff = 0.25; % Hz
num_poles = 4;

[lowpassB, lowpassA] = butter(num_poles, low_pass_cutoff/(sps/2), 'low');

impulse_response = impz(lowpassB, lowpassA);

plot(impulse_response);
title('Impulse Response of Low-pass Butterworth Filter');
xlabel('Sample');
ylabel('Amplitude');

%% 





%% 
% Compute power spectra
[pxx_original,f_original] = periodogram(HzWindSpd_T,[],[],sps);
[pxx_lowp,f_lp] = periodogram(wind_speed_lowpass,[],[],sps);
[pxx_highp,f_hp] = periodogram(wind_speed_highpass,[],[],sps);
[pxx_bandp,f_bp] = periodogram(wind_speed_bandpass,[],[],sps);

% Apply Hann window for spectral smoothing
window_length = 21; % window length
window = hann(window_length);
window = window/sum(window);

pxx_orig_smooth = conv(pxx_original,window,'same');
pxx_lp_smooth = conv(pxx_lowp,window,'same');
pxx_hp_smooth = conv(pxx_highp,window,'same');
pxx_bp_smooth = conv(pxx_bandp,window,'same');

% Plot power spectra

semilogx(f_original,10*log10(pxx_orig_smooth),'k','LineWidth',1)
hold on

semilogx(f_lp,10*log10(pxx_lp_smooth),'b','LineWidth',1)
semilogx(f_hp,10*log10(pxx_hp_smooth),'r','LineWidth',1)
semilogx(f_bp,10*log10(pxx_bp_smooth),'m','LineWidth',1)

ylabel('Power (dB/Hz)')
xlabel('Frequency (Hz)')
title('Power Spectra of Wind Speed from Lower Deer Point, Idaho')
legend('Original Data','Low-pass Filtered','High-pass Filtered','Band-pass Filtered')
grid on

%%
% Compute power spectra
[pxx_original,f_original] = periodogram(HzWindSpd_T,[],[],sps);
[pxx_lowp,f_lp] = periodogram(wind_speed_lowpass,[],[],sps);
[pxx_highp,f_hp] = periodogram(wind_speed_highpass,[],[],sps);
[pxx_bandp,f_bp] = periodogram(wind_speed_bandpass,[],[],sps);

% Apply Hann window for spectral smoothing
window_length = 21; % window length
window = hann(window_length);
window = window/sum(window);

pxx_orig_smooth = conv(pxx_original,window,'same');
pxx_lp_smooth = conv(pxx_lowp,window,'same');
pxx_hp_smooth = conv(pxx_highp,window,'same');
pxx_bp_smooth = conv(pxx_bandp,window,'same');

% Plot power spectra
figure;
semilogx(f_original,10*log10(pxx_orig_smooth),'k','LineWidth',1)
hold on
semilogx(f_lp,10*log10(pxx_lp_smooth),'b','LineWidth',1)
% semilogx(f_hp,10*log10(pxx_hp_smooth),'r','LineWidth',1)
% semilogx(f_bp,10*log10(pxx_bp_smooth),'m','LineWidth',1)
% ylabel('Power (dB/Hz)')
xlabel('Frequency (Hz)')
title('Power Spectra of Wind Speed from Lower Deer Point, Idaho')
legend('Original Data','Low-pass Filtered','High-pass Filtered','Band-pass Filtered')
grid on



%%

sps = 2000; 
ts = (1:length(HzWindSpd_T))/sps;

vel = detrend(HzWindSpd_T, 'constant');

% create spectrogram
window_size = 1024; % choose window size
overlap_ratio = 0.5; % choose overlap ratio
nfft = window_size; % use window size for FFT length
[s,f,t] = spectrogram(vel,window_size,round(window_size*overlap_ratio),nfft,sps);

% plot spectrogram
imagesc(t,f,10*log10(abs(s)));
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('Spectrogram of Wind Speed from Lower Deer Point, Idaho');
colorbar;

hold on;
plot(ts,vel,'k');
colormap(jet)

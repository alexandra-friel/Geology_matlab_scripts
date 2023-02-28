path = 'C:/Users/lexiefriel/Desktop'; 
load whistle.m

downsampling = sound_variable(1:12:end);  
plot(downsampling);

xlabel('Sample');
ylabel('Amplitude');
nyquist =  frequency/(2*12); %12 from 1:12:end

sound(downsampled, frequency/12);

%% 
resampled_sound = decimate(sound_variable,12); %I need to download the signal processing tool box. 

%% Short answer

disp('The whistles will vary in sound because it was measured at a value lower than the Nyquist rate.'); 
disp('This is is the minimum sample rate that is necessary to avoid potential aliasing. ')

disp('The decimate function will reduce the sample in order to remove high frequency components that may attribute to aliasing.')
disp('This may be a useful anti-aliasing filter. ')






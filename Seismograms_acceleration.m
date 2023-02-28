%% Name: Alexandra Friel
%Date: February 15th, 2023
%Digital Signal Processing Homework 03
%This script is designed to read in a data set from Mt. Erebus, a volcano
%from Antartica. It will plot the acceleration of the seismogram and the
%scaled displacement. 


%Variables given: data, atod, hdr, sps, hdr.atod, hdr.sps

%%
directory = 'C/Users/lexiefriel/Desktop'; 
load('Erebus_seismogram.mat') %Given .mat file. 

time = (0:length(data)-1)/hdr.sps;
sensitivity = 3200;
atod = hdr.atod;

data = detrend(data);  %will detrend the given data.
data= (data * atod / sensitivity) * 1000; %will not plot if it is not called data. 

plot(time, data,'m')
xlabel('Time (seconds)')
ylabel('Amplitude (mm/sec)')
title('Seismogram from Mount Erebus')

%% Properly scaled acceleration seismogram

dt = 1/hdr.sps; 
acceleration = diff(data)/dt * sensitivity * 1000; 
data = detrend(data); 

time = (0:numel(acceleration)-1) * dt; %Numel -> number of elements in an array/expression
plot(time, acceleration,'r'); % Plot the acceleration seismogram
xlabel('Time (s)');
ylabel('Acceleration (mm/s^2)');
title('Seismic activity from Mount Erebus')

%% Properly scaled displacement seismogram

velocity = data * sensitivity * 1000; 
displacement = zeros(size(velocity));
displacement(2:end) = cumsum(0.5 * (velocity(2:end) + velocity(1:end-1)) * dt * 1000);
time = (0:numel(acceleration)-1) * dt; 

plot(t, displacement,'r'); 
xlabel('Time (s)');
ylabel('Displacement (mm)');
title('Seismic activity from Mount Erebus')


%% Subplot

subplot(2,2,1)
plot(t, displacement,'m'); 
xlabel('Time (s)');
ylabel('Displacement (mm)');
title('Seismic activity from Mount Erebus')

subplot(2,2,2)
plot(time, acceleration,'r'); 
xlabel('Time (s)');
ylabel('Acceleration (mm/s^2)');
title('Seismic activity Mount Erebus')

subplot(2,2,3)
plot(t, displacement,'k'); 
xlabel('Time (s)');
ylabel('Displacement (mm)');

title('Seismic activity from Mt. Erebus')

%title_handle = title('Visualizing seismic activity from Mt. Erebus'); 
%set(title_handle, 'Position', [-2,4.100,10])

%set(gcf, 'Position', [100 100 1000 800]);



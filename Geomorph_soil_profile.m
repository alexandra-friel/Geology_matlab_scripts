%% This script was written by Alexandra Friel for Geomorphology Soil Analysis 

%%

depth = [4 10 54 140 150 180]; %O layer to C layer, in cm

sand = [15 20 25 45 60 65]; %percent sand
silt = [1 10 11 25 30 32]; %percent silt
clay = [80 70 64 30 10 3]; %percent clay

depth_smooth = linspace(depth(1), depth(end), 1000);
sand_smooth = spline(depth, sand, depth_smooth);
silt_smooth = spline(depth, silt, depth_smooth);
clay_smooth = spline(depth, clay, depth_smooth);

figure;
plot(depth_smooth, sand_smooth, '-.k', 'LineWidth', 2);
hold on;
plot(depth_smooth, silt_smooth, '-.m', 'LineWidth', 2);
plot(depth_smooth, clay_smooth, '-.b', 'LineWidth', 2);
ylabel('Percentage (%)');
xlabel('Depth (cm)');
title('Sand, silt, & clay percentages for Gowen Terrace Soil Profile')
legend('Sand', 'Silt', 'Clay');

% Set axis limits
xlim([0 180]);
ylim([0 100]);


%% Percentage of Calcium Carbonate

depth = [4, 10, 54, 140, 150, 180]; % depth in cm
CaCO3_percent = [0, 0, 5, 10, 60, 80]; % CaCO3 percentage

CaCO3_smooth = movmean(CaCO3_percent, window_size);

plot(depth, CaCO3_smooth,'m'); % plot smoothed curve
hold off;

% Add labels and title
xlabel('Depth (cm)');
ylabel('CaCO3 percentage (%)');
title('CaCO3 percentages vs. Depth at Gowen Terrace');

%% Gravel 
depth = [4 10 54 140 150 180];
gravel_percentage = [1 5 10 15 25 30]; 

plot(depth, gravel_percentage,'b')
xlabel('Depth in cm')
ylabel('Gravel content percentage')
title('percentage of gravel vs. depth')

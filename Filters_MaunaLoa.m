cd('/Users/lexiefriel/Desktop');
load('CO2_MaunaLoa.mat')

%%
CO2_filtered = medfilt1(CO2,3); 
figure(1)
plot(ts,CO2_filtered, "o",'MarkerSize',3,'Color','magenta')
hold on; 

plot(ts,CO2_filtered,"-",'MarkerSize',2,'Color','blue')

title('Concentration of CO2 vs. time using Median Filter')
xlabel('time in years')
ylabel('Concentration of CO2 in ppm')

grid on; 
hold off; 

%%
n_points = 31;
theta_range = [0, pi];
theta = linspace(theta_range(1), theta_range(2), n_points);

sin_theta = sin(theta);
sin_theta_norm = sin_theta / sum(sin_theta);
filtered_CO2 = conv(CO2_filtered, sin_theta_norm, 'valid');
%% 

figure;
plot(ts, CO2_filtered, 'k');
hold on;
plot(ts(16:end-15), filtered_CO2, 'm');

hold off; 
title('Concentration of CO2 vs. time using convolution filter')
xlabel('time in years')
ylabel('Concentration of CO2 in ppm')
legend('Original signal', 'Filtered signal');
grid on;

%% 

filt_fig6_3b = -sin(0:pi/30:pi);
filt_fig6_3b = filt_fig6_3b / sum(filt_fig6_3b) * -1;
% theta = 0:pi/30:pi;
% filt_fig6_3b = -sin(theta);
% filt_fig6_3b = filt_fig6_3b / sum(abs(filt_fig6_3b)) * -1;

figure(3)
plot(ts,CO2_filtered)
hold on; 
filt_fig6_3b = sum(filt_fig6_3b(1:15)* -2); 

output2 = conv(CO2_filtered, filt_fig6_3b); 
filtered_outp2 = output2(30:length(CO2_filtered)); 

plot(ts(15:end-15),filtered_outp2,'m')
title('Concentration of CO2 vs. time using ')
xlabel('time in years')
ylabel('Concentration of CO2 in ppm')
hold off; 
legend('Original signal', 'Filtered signal');
grid on;


%% 
testing = conv(filt_fig6_3a, filt_fig6_3b); 
figure(4) 
plot(testing,'-.m')
title('Combined filters from A & B'); 
legend('A and B')
grid on;

%%
figure(5)
combined_data = filtered_outp + filtered_outp2; 
plot(ts,CO2_filtered)

hold on; 
start_index = 15;
end_index = length(ts) - 15;
plot(ts(start_index:end_index), combined_data,'m')
title('Normal vs. high pass filtered trend')
legend('Normal', 'high pass filtered signal');
hold off; 


%%
%The high pass filter would remove the DC offset. 
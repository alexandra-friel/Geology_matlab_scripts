% Import data from Excel sheet
data = readtable('pebblecount.xlsx');

%%
pebbleData = data{:,1};
mu = mean(pebbleData);
sigma = std(pebbleData);
x = mu - 4*sigma : 0.01 : mu + 4*sigma;
y = normpdf(x, mu, sigma);
plot(x, y, 'Color', [161 135 245]/255, 'LineWidth',3);
xlabel('Pebble Count');
ylabel('Probability Density');
title('Pebble Count Bell Curve');

%%
figure('Color',[233 226 253  ]/255);
histogram(pebbleData, 'FaceColor', [161 135 245]/255, 'EdgeColor', [104, 88, 156]/255, 'LineWidth',2);
title('Pebble counts at CottonWood Creek')
xlabel('Pebble B-axis measurement')
%%
meanPebble = mean(pebbleData); 
medianPebble = median(pebbleData); 
modePebble = mode(pebbleData); 
rangePebble = range(pebbleData); 

fprintf('Mean Pebble Data: %.2f\n', meanPebble);
fprintf('Median Pebble Data: %.2f\n', medianPebble);
fprintf('Mode Pebble Data: %.2f\n', modePebble);
fprintf('Range Pebble Data: %.2f\n', rangePebble);
fprintf('\n')


%% Measure assymetry of distributiono
skewnessPebble = skewness(pebbleData);
fprintf('Skewness of Pebble Data: %.2f\n', skewnessPebble);
fprintf('\n')




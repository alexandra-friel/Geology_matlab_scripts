%% Frequency domain manipulation - cross correlation

f = [0 0 1 -1 0 0 0]; 
g = [0 0 0 2 -2 0 0]; 

co = conv(f,g); 
xc = xcorr(f,g); 

%Convert to freq domain

F = fft(f);
G = fft(g);

figure(1)
subplot(3,1,1)
stem(co,'Filled','m')
title('Time domain convolution')

subplot(3,1,2)
CO = F.*G;  % . for two vectors multiplied together. It is now a complex vector
stem(ifft(CO),'filled','b') %one piece of a repeating time series. Very periodic time series. Sinusoids will go onto infinity. Takes complex numbers and brings them back into the time domain.
title('') 

subplot(3,1,3)
CO = F .*G; 
stem(fftshift(ifft(CO)),'Filled','r')

%'same' will truncate zeros 

%%
figure(2); clf;
xc = xcorr(f,g); 

subplot(2,1,1)
stem(xc,'m','filled')

subplot(2,1,2)
XC = F.*conj(G); %Freq domain
stem(-3:3,fftshift(ifft(XC)),'b','filled')



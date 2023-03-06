%% DSP in-class work: March 6th, 2023

t = (0:63); 
x = t*0; 
x(4) = 1; 

subplot(3,1,1)
plot(t,x,'m')

subplot(3,1,2)
plot(t,abs(fft(x)),'k')

subplot(3,1,3)
plot(t,angle(fft(x)),'b')


%% Frequency domain representations 

figure(2) 

t = (-20:.1:20); 
x = sinc(t);  % sinc = sin(pi*x) / (pi)(x)

subplot(2,1,1)
plot(t,x,'om','MarkerSize',2)

X = fft(x); 
subplot(2,1,2)
plot(abs(X),'-r')

%% Gaussian pulse - normal distribution 

figure(3); 
t = (-20:.1:20); 
sigma = 1; 
x = 1/sigma/sqrt(2*pi)*exp(-.5*((t/sigma).^2)); 

subplot(2,1,1)
plot(t,x,'m')

subplot(2,1,2) %convert to frequency domain
X = fft(x) 
plot(abs(X),'k') %Centered at 0 frequency


figure(4); 
t = (-20:.1:20); 
x = t*0; 
x(10:100:end)=1; 
subplot(2,1,1)
plot(t,x,'m')

subplot(2,1,2)
plot(abs(fft(x)),'m')



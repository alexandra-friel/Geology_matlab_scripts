npoles = 2;
corner = 1; 
sps = 100; 

[B,A] = butter(npoles, [corner]/(sps/2),'low')

%Recursive nature allows minimal coefficients for infinite impulse
%response. 
%Finite: h[n] = [1 -1]
% Can do this with only a few values 

%% 

t0 = 6.655 %time constant 
fc = 1/2/pi/t0; 
sps = 100; 
corner = fc*sps %Retain /throw out information
npoles = 1; 

[B,A] = butter(npoles, [corner]/(sps/2),'low') %Low pass filtering


%% 

% corner (hz)  = fc * sample rate 



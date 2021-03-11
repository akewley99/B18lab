%Section B EEG
%% Exercise 1
%question 1
%load data
EEG_closed=load('C:\Users\alexk\OneDrive\Documents\MATLAB\B18-Wearables-Laboratory-main\B18-Wearables-Laboratory-main\DATA\B18_EEG_data\EEGeyesclosed.mat');
EEG_open=load('C:\Users\alexk\OneDrive\Documents\MATLAB\B18-Wearables-Laboratory-main\B18-Wearables-Laboratory-main\DATA\B18_EEG_data\EEGeyesopen.mat');

closedfull=EEG_closed.eyesclosed;
openfull=EEG_open.eyesopen;

closed=closedfull(:,1:512);
open=openfull(:,1:512);


%convert samples to time
fs =256; % Sampling frequency (Hz)

samples_c=1:length(closed); % 'signal' is your vector of ecg values
%'time' is a corresponding monotonic vector of time values
time = samples_c./fs;

%plot
figure
subplot(2,1,1)
plot(time,closed)
title('Eyes Closed')
xlabel('Time/s')
ylabel('Amplitude/\muV')
%ylim([-200 200])

subplot(2,1,2)
plot(time,open)
title('Eyes Open')
xlabel('Time/s')
ylabel('Amplitude/\muV')
%ylim([-200 200])


%question 2 
%low pass filter, detrend and normalise

%detrend by subtracting mean
D_closed=closed-mean(closed);
D_open=open-mean(open);

%standardise using z-score
N=length(closed);
sigma_c=(sum(D_closed.^2)/(N-1))^0.5;
sigma_o=(sum(D_open.^2)/(N-1))^0.5;

Z_closed=D_closed./sigma_c;
Z_open=D_open./sigma_o;

%apply low pass filter
%determine parameters by reporting and 
%analysing various combos

fc=35;       %cut off freq in Hz must be <256
norder=4;   %filter order
Wn=fc/fs;   %cut off freq

%design filter
[B,A]=butter(norder,Wn);

%apply LPF
filtered_closed=filter(B,A,Z_closed);
filtered_open=filter(B,A,Z_open);

%plot filtered samples
figure
subplot(2,1,1)
plot(time,filtered_closed)
ylabel('Amplitude/\muV')
xlabel('Time/s')
title('Closed After Filtering')

subplot(2,1,2)
plot(time,filtered_open)
ylabel('Amplitude/\muV')
xlabel('Time/s')
title('Open After Filtering')


%plot modified periodogram up to 60Hz
%hamming window function
window_c=hamming(length(filtered_closed));
window_o=hamming(length(filtered_open));

%modified periodogram
[pxx_c,f_c]=periodogram(filtered_closed,window_c,[],fs,'psd');
[pxx_o,f_o]=periodogram(filtered_open,window_o,[],fs,'psd');


%plot periodogram
figure
subplot(1,2,1)
plot(f_c,pxx_c);
xlim([0,60])
xlabel('Frequency/Hz')
ylabel('Power Spectral Density/mV^{2}Hz^{-1}')
title('Modified Periodogram (closed)')
%ylim([0,5])

subplot(1,2,2)
plot(f_o,pxx_o);
xlim([0,60])
xlabel('Frequency/Hz')
ylabel('Power Spectral Density/mV^2Hz^-1')
title('Modified Periodogram (open)')
%ylim([0,5])





%% Exercise 2
%question 1
%load sleep state data
EEG_SSD=load('C:\Users\alexk\OneDrive\Documents\MATLAB\B18-Wearables-Laboratory-main\B18-Wearables-Laboratory-main\DATA\B18_EEG_data\EEGSleepStateData.mat');

SSD=EEG_SSD.EEGSleepStateData; 

%question 2
%plot hypnogram
num_cat=SSD{:,3};           %numerical category
timeSSD=1:length(num_cat);  %time vector

%plot
figure
stairs(timeSSD,num_cat)
xlabel('Time/mins')
ylabel('Sleep State')
ylim([0.5 6.5])
title('Hypnogram of Sleep Categories')


%question 3 
%plot periodogram PSD
%detrend
D_num_cat=num_cat-mean(num_cat);
fs=60;  %sampling rate

%obtain periodogram
window=hamming(length(D_num_cat));
[pxx,f]=periodogram(D_num_cat,window,[],fs);

%plot
figure
plot(f,pxx);
xlim([0,30])
xlabel('Frequency/cycles per hour')
ylabel('Power Spectral Density/mV^{2}Hz^{-1}')
title('Modified Periodogram SSD')







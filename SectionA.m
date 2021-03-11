%Section A ECG
%% Excercise 1
%load ECG data
ECG=load('C:\Users\alexk\OneDrive\Documents\MATLAB\B18-Wearables-Laboratory-main\B18-Wearables-Laboratory-main\DATA\B18_ECG_data\PhysionetData.mat');


%question 1
%load data into tables
RecordingInfo=ECG.RecordingInfo; %struct can be addressed using . 
Signals=ECG.Signals;
tabulate(RecordingInfo.Label)


%question 2
%select signals
normalECG=Signals{7}; %cell needs curly {}
normalECG=normalECG/100;

%convert samples to time
fs =300; % Sampling frequency (Hz)
samples=1:length(normalECG); % 'signal' is your vector of ecg values
%'time' is a corresponding monotonic vector of time values
time = samples./fs;

%plot normal ECG
figure
X1=time(:,100:500);
Y=normalECG(:,100:500);
plot(X1,Y)
ylabel('Amplitude/mV')
xlabel('Time/s')
text(0.53,1.25,'P')
text(0.61,-0.15,'Q')
text(0.59,3.80,'R')
text(0.71,-1.55,'S')
text(0.86,1.14,'T')
print -depsc epssecAex1q2



%question 3

%select signals
%already have normal
afECG=Signals{6}/100;
otherECG=Signals{15}/100;
noisyECG=Signals{197}/100;

%select points
X2=time(:,1:6000);
NOR=normalECG(:,1:6000);
AF=afECG(:,1:6000);
O=otherECG(:,1:6000);
NOI=noisyECG(:,1:6000);

%plot
figure
subplot(4,1,1)
plot(X2,NOR)
ylabel('Amplitude/mV')
xlabel('Time/s')
title('Normal Rhythm; Recording: A01021')

subplot(4,1,2)
plot(X2,AF)
ylabel('Amplitude/mV')
xlabel('Time/s')
title('AF Rhythm; Recording: A08003')

subplot(4,1,3)
plot(X2,O)
ylabel('Amplitude/mV')
xlabel('Time/s')
title('Other Rhythm; Recording: A08511')

subplot(4,1,4)
plot(X2,NOI)
ylabel('Amplitude/mV')
xlabel('Time/s')
title('Noisy; Recording: A02514')

print -depsc epssecAex1q3



%% Exercise 2
%use normalECG and afECG from earlier

%question 1

%get time again
%convert samples to time
fs =300; % Sampling frequency (Hz)
samples=1:length(normalECG); % 'signal' is your vector of ecg values
%'time' is a corresponding monotonic vector of time values
time = samples./fs;

%detrend by subtracting mean
DnormalECG=normalECG-mean(normalECG);
DafECG=afECG-mean(afECG);

%standardise using z-score

%ZnormalECG=(DnormalECG./std(normalECG))';
%ZafECG=(DafECG./std(afECG))';

N=length(normalECG);
sigma_n=(sum(DnormalECG.^2)/(N-1))^0.5;
sigma_af=(sum(DafECG.^2)/(N-1))^0.5;

ZnormalECG=DnormalECG./sigma_n;
ZafECG=DafECG./sigma_af;


%question 2
%apply band pass filtering
fs=300;
f_low=40;
f_high=5;
norder=4;

filtered_normal=filter_ecg(ZnormalECG,fs,f_high,f_low,norder);
filtered_af=filter_ecg(ZafECG,fs,f_high,f_low,norder);




%question 3 
%plotq1&2

figure
subplot(6,1,1)
plot(time(:,2000:6000),DnormalECG(:,2000:6000))
ylabel('Amplitude/mV')
xlabel('Time/s')
title('Normal After Detrending')

subplot(6,1,4)
plot(time(:,2000:6000),DafECG(:,2000:6000))
ylabel('Amplitude/mV')
xlabel('Time/s')
title('AF After Detrending')

subplot(6,1,2)
plot(time(:,2000:6000),ZnormalECG(:,2000:6000))
ylabel('Amplitude/mV')
xlabel('Time/s')
title('Normal After Standardisation')

subplot(6,1,5)
plot(time(:,2000:6000),ZafECG(:,2000:6000))
ylabel('Amplitude/mV')
xlabel('Time/s')
title('AF After Standardisation')

subplot(6,1,3)
plot(time(:,2000:6000),filtered_normal(:,2000:6000))
ylabel('Amplitude/mV')
xlabel('Time/s')
title('Normal After Filtering')

subplot(6,1,6)
plot(time(:,2000:6000),filtered_af(:,2000:6000))
ylabel('Amplitude/mV')
xlabel('Time/s')
title('AF After Filtering')


print -depsc epssecAex2q3


%% Exercise 3
%question 1
%plot periodogram PSD
%use filtered_normal and filtered_af

%hamming window function
window_n=hamming(length(filtered_normal));
window_af=hamming(length(filtered_af));

%modified periodogram
[pxx_n,f_n]=periodogram(filtered_normal,window_n,[],fs,'psd');
[pxx_af,f_af]=periodogram(filtered_af,window_af,[],fs,'psd');


%plot periodogram
figure
subplot(1,2,1)
plot(f_n,pxx_n);
xlim([0,f_low])
xlabel('Frequency/Hz')
ylabel('Power Spectral Density/mV^{2}Hz^{-1}')
title('Modified Periodogram (normal)')
%ylim([0,5])

subplot(1,2,2)
plot(f_af,pxx_af);
xlim([0,f_low])
xlabel('Frequency/Hz')
ylabel('Power Spectral Density/mV^2Hz^-1')
title('Modified Periodogram (AF)')
%ylim([0,5])

print -depsc epssecAex3q1

%welch's PSD
[pxx_pwelch_n,f_pwelch_n]=pwelch(filtered_normal,[],[],[],fs,'psd');
[pxx_pwelch_af,f_pwelch_af]=pwelch(filtered_af,[],[],[],fs,'psd');

%plot welchs
figure
subplot(1,2,1)
plot(f_pwelch_n,pxx_pwelch_n);
xlim([0,f_low])
xlabel('Frequency/Hz')
ylabel('Power Spectral Density/mV^{2}Hz^{-1}')
title('Welchs PSD (normal)')

subplot(1,2,2)
plot(f_pwelch_af,pxx_pwelch_af);
xlim([0,f_low])
xlabel('Frequency/Hz')
ylabel('Power Spectral Density/mV^2Hz^-1')
title('Welchs PSD (AF)')


%% Exercise 4
%Section 1
%using filtered_normal & filtered_af
%question 1
%r peak with findpeaks.m

%time=X2;
MinPeakHeight=1.4;
MinPeakDistance=2;

[R_pks_n,R_locs_n]=findpeaks(filtered_normal,...
    'MinPeakHeight',MinPeakHeight,...
    'MinPeakDistance',MinPeakDistance);
[R_pks_af,R_locs_af]=findpeaks(filtered_af,...
    'MinPeakHeight',MinPeakHeight,...
    'MinPeakDistance',MinPeakDistance);

figure
subplot(2,1,1)
plot(time,filtered_normal); hold on;
plot(time(R_locs_n),filtered_normal(R_locs_n),'ro')
xlabel('Time/s')
ylabel('Amplitude/mV')
title('R Peak Detection (Normal)')

subplot(2,1,2)
plot(time,filtered_af); hold on;
plot(time(R_locs_af),filtered_af(R_locs_af),'ro')
xlabel('Time/s')
ylabel('Amplitude/mV')
title('R Peak Detection (AF)')

%question 2
%determine mean heart rate
%=r peaks per time

[~,a]=size(R_pks_n);
[~,b]=size(R_pks_af);

time_R_n=time(R_locs_n);
time_R_af=time(R_locs_af);

[~,c]=size(time_R_n);
end_time_n=time_R_n(1,c);

[~,d]=size(time_R_af);
end_time_af=time_R_af(1,d);

heartrate_n=a/end_time_n;
heartrate_af=b/end_time_af;

%in bpm
bpm_n=heartrate_n*60;
bpm_af=heartrate_af*60;

display(bpm_n)
display(bpm_af)




%Section 2

%question 3
%get R-R interval times

RR_n=1000*(diff(R_locs_n)./fs);
RR_af=1000*(diff(R_locs_af)./fs);

%plot RR intervals over time
figure
subplot(2,1,1)
plot(time_R_n(:,2:c),RR_n,'--o')
xlabel('Time/s')
ylabel('R-R Interval/ms')
title('R-R Interval Times (Normal)')
ylim([400 1000])

subplot(2,1,2)
plot(time_R_af(:,2:d),RR_af,'--o')
xlabel('Time/s')
ylabel('R-R Interval/ms')
title('R-R Interval Times (AF)')
ylim([400 1000])


%question 4
%plot histogram
Binwidth=20;
figure
subplot(2,1,1)
histogram(RR_n,'Binwidth',Binwidth)
xlabel('R-R Interval/ms')
ylabel('Count')
xlim([400,1100])
ylim([0,9.5])
title('Distribution of R-R Intervals (Normal)')

subplot(2,1,2)
histogram(RR_af,'Binwidth',Binwidth)
xlabel('R-R Interval/ms')
ylabel('Count')
xlim([400,1100])
ylim([0,9.5])
title('Distribution of R-R Intervals (AF)')


%question 5
%plot lomb_scargle periodogram
%RR interval in seconds
RRi_n=RR_n/1000;
RRi_af=RR_af/1000;

%detrend
RRi_n=RRi_n-mean(RRi_n);
RRi_af=RRi_af-mean(RRi_af);

%time from first interval
RRit_n=time(R_locs_n(2:end));
RRit_af=time(R_locs_af(2:end));

%frequencies to be evaluated
f_interest=0.001:0.001:0.6;

%compute lomb-scargle periodogram
[pxx_plomb_n,f_plomb_n]=plomb(RRi_n,RRit_n,f_interest,'psd');
[pxx_plomb_af,f_plomb_af]=plomb(RRi_af,RRit_af,f_interest,'psd');

%plot PSD
figure
subplot(2,1,1)
plot(f_plomb_n,pxx_plomb_n)
xlabel('Frequency/Hz')
ylabel('Power Spectral Density/s^{2}Hz^{-1})')
title('Lomb-Scargle PSD Spectrum (Normal)')

subplot(2,1,2)
plot(f_plomb_af,pxx_plomb_af)
xlabel('Frequency/Hz')
ylabel('Power Spectral Density/s^{2}Hz^{-1})')
title('Lomb-Scargle PSD Spectrum (AF)')

%question 6
%find HRV measures
HRV_n=hrv_measures(RR_n,f_plomb_n,pxx_plomb_n);
HRV_af=hrv_measures(RR_af,f_plomb_af,pxx_plomb_af);

display(HRV_n)
display(HRV_af)









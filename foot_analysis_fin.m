% foot_analysis_fin.m
% First peak can be ignored as it is due to buffer function padding zeroes
% in the beginning.

%% Loading Dataset
clc
clear all
close all
[s, h] = sload('gautam-real-2-[2016.01.21-13.04.18].gdf'); %Load GDF File
fs = h.SampleRate;

%% Select the Channel and Band for Filtering
x = find(ismember(h.Label,['Cz']));
signal = (s(:,x));                 %Cz
passBand.low = 12;
passBand.high = 30;

order = 4; 
lowFreq = passBand.low * (2/fs);
highFreq = passBand.high * (2/fs);
[B A] = butter(order, [lowFreq highFreq]);
signal = filter(B, A, signal);

%% Plotting the Beta-Filtered Signal and Windowing
t = 0 : 1/fs : (length(signal) - 1) / fs;
figure;
plot(t,signal); xlabel('Time(s)'); ylabel('Amplitude(uV)'); title('Beta-Filtered Signal(16Hz-24Hz)');
axis([0 150 -50 50]);

overlap_factor = 0.9;
mat = buffer(signal, fs, ceil(overlap_factor * fs));           %create a matrix with overlapped segments, also called 'Time Based Epoching'. 1 sec segments every 0.1 secs.
mat = mat';                                                    %segments created by buffer are column vectors, so take the transpose.
mat = mat .* mat;                                              %square the signal. (Simple DSP = x * x)

signal_avg = mean(mat, 2);                                     %mean of the segments


window_marker = ceil(h.EVENT.POS / ((1 - overlap_factor) * fs));

%% Computing Epoch Average and Dataset Statistics
for i = 1:size(signal_avg,1) - 3
    epoch_avg(i) = mean(signal_avg(i:i+3));
end

calib_segment = epoch_avg(window_marker(1) : window_marker(2));
calib_mean = mean(calib_segment)
calib_var = var(calib_segment)
calib_sd = sqrt(calib_var)

%% Creation of Event Markers
j = 1;
for i = 1:length(epoch_avg)
    event_marker(i) = 0;
    if(i == window_marker(j)&& j<length(window_marker))
        event_marker(i) = 5*calib_mean;
        j = j + 1;
    end
end

new = mean(epoch_avg);
th_m(1 : length(epoch_avg), 1) = mean(epoch_avg);

figure;
subplot(2, 1, 1);
plot(epoch_avg); title('Epoch-Average-Power');
axis([0 1500 -10 5*calib_mean]);
hold on
stem(event_marker);
subplot(2,1,2);
spectrogram(signal, fs, overlap_factor*fs, fs, fs,'yaxis');
colormap hsv

%% Calib_min and Calib_max need to be manually configured on a per signal
%basis.
calib_min = calib_mean + 3*calib_sd;
calib_max = calib_mean + 6*calib_sd;


%% Cropping the signal between (Mean + 3*Variance) and (Mean + 6*Variance)
for i = 1:size(epoch_avg,2)
    if(epoch_avg(i) > calib_max)                      
       epoch_avg(i) = calib_max;
    end
    if(epoch_avg(i) < calib_min)
        epoch_avg(i) = calib_min;
    end
end

%% Center around (Mean + 3*Variance)
beta_power = epoch_avg - (calib_max + calib_min)/2;

th(1:length(beta_power),1) = calib_min;   
th(1:length(beta_power),1) = 0;

%% Plotting the Beta-Power Cropped Signal
figure;
plot(beta_power); title('Beta-Power-Cropped with Threshold');
%axis([0 1500 -10 20]);
hold on
plot(th','r');

%% Classifying Based on Threshold
class_true = 0;
class_false = 0;
for cl = 1:length(beta_power) - 1
	if ((beta_power(cl) < th(1,1)) && (beta_power(cl+1) > th(1,1)))
		class_true = class_true + 1;
end
end

disp('The no. of true classes: ');
disp(class_true-1);                                     %first classification to be ignored
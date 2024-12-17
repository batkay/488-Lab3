load("NeuroDOT_Data_Sample_OUT1.mat");

figure;
t = 1:size(data, 2);
y = log(bsxfun(@times,data,1./mean(data,2)));


% out = y(1:size(y, 1)/2, :);
out = y;
out = out((info.pairs.NN==2 | info.pairs.NN == 1) & info.pairs.WL == 2, :);


subplot(2, 2, 1);
plot(t, out);
hold on;
for i = 1:length(info.paradigm.Pulse_2), xline(info.paradigm.synchpts(info.paradigm.Pulse_2(i))), end



deviation = std(out, 0, 2);
threshold = 0.1;
valid_rows = deviation <= threshold; % Identify rows with std <= threshold
data_wave_1 = out(valid_rows, :); % Retain only valid rows

subplot(2, 2, 2);
plot(t, data_wave_1);
hold on;
for i = 1:length(info.paradigm.Pulse_2), xline(info.paradigm.synchpts(info.paradigm.Pulse_2(i))), end



% valid_rows = max(abs(data_wave_1), [], 2) <= threshold; % Identify rows with std <= threshold
% data_wave_1 = data_wave_1(valid_rows, :); % Retain only valid rows

subplot(2, 2, 3);
plot(t, data_wave_1);
hold on;
for i = 1:length(info.paradigm.Pulse_2), xline(info.paradigm.synchpts(info.paradigm.Pulse_2(i))), end



% [filtered, D] = lowpass(data_wave_1', 0.003, 11, StopbandAttenuation=60, Steepness=.5, ImpulseResponse="fir");
% [filtered, D] = bandpass(data_wave_1', [0.00000003, 1], 11, StopbandAttenuation=60);
% [filtered] = lowpass(data_wave_1', 0.1, 11, StopbandAttenuation=60, Steepness=.5, ImpulseResponse="fir");
% [filtered, D] = highpass(data_wave_1', 1, 11);
[B, BatA] = butter(2, [.02, 1]/(11/2), "bandpass");
filtered = filter(B, BatA, data_wave_1);
filtered = movmean(filtered', 10)';
% filtered = filtered';

subplot(2, 2, 4);
plot(t, filtered);
hold on;
for i = 1:length(info.paradigm.Pulse_2), xline(info.paradigm.synchpts(info.paradigm.Pulse_2(i))), end


fs = 11;
ff = [0:fs/size(data_wave_1,2):fs/2]; % frequency vector 
Y = fft(data_wave_1, [], 2);
Yf = fft(filtered,[],2); % fft along 2nd dimension of y, all measurements 


figure;
subplot(1, 2, 1);
plot(ff,abs(Y(:, 1:length(ff))));
subplot(1, 2, 2);
plot(ff,abs(Yf(:, 1:length(ff))));

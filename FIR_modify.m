close all;
clear;
clc;
data = readtable('C:/Users/Mega Store/Desktop/100_signals.csv');
ecg = data.MLII;
Fs = 360;
N = 3600;
ecg_segment = ecg(1:N);
t_segment = (0:N-1) / Fs;
ord_fir = 500;

% High pass
b_hp = fir1(ord_fir,0.5/(Fs/2),'high');
a_hp = 1;

% FIR Notch Filter
b_notch = fir1(ord_fir,[49 50]/(Fs/2),'stop');
a_notch = 1;

% Low Pass
b_lp = fir1(ord_fir,100/(Fs/2),'low');
a_lp = 1;


ecg_hp = filtfilt(b_hp, a_hp, ecg_segment);
ecg_notch = filtfilt(b_notch, a_notch, ecg_hp);
ecg_clean = filtfilt(b_lp, a_lp, ecg_notch);

figure;
[pxx_orig,  f] = pwelch(ecg_segment, [], [], [], Fs);
[pxx_clean, ~] = pwelch(ecg_clean,   [], [], [], Fs);
plot(f, 10*log10(pxx_orig),'r', 'LineWidth', 1.2);
hold on;
plot(f, 10*log10(pxx_clean),'b', 'LineWidth', 1.2);
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density (dB/Hz)');
title('PSD: Original vs Filtered ECG');
legend('Original', 'Filtered');
xlim([0 180]);
grid on;

figure;
subplot(2,1,1);
spectrogram(ecg_segment, 256, 200, 256, Fs, 'yaxis');
title('Original ECG Spectrogram');

subplot(2,1,2);
spectrogram(ecg_clean, 256, 200, 256, Fs, 'yaxis');
title('Filtered ECG Spectrogram');

figure;
subplot(2,1,1);
plot(t_segment, ecg_segment,'r');
title('Original ECG');
ylabel('Amplitude (mV)');
grid on;

subplot(2,1,2);
plot(t_segment, ecg_clean,'b');
title('Filtered ECG');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
grid on;

% SNR Report
idx_baseline  = f < 0.5;
idx_powerline = f >= 49  & f <= 51;
idx_emg = f >= 100 & f <= 150;
idx_signal = f >= 0.5 & f <= 100;
noise_power_before = sum(pxx_orig(idx_baseline)) + sum(pxx_orig(idx_powerline)) + sum(pxx_orig(idx_emg));
noise_power_after  = sum(pxx_clean(idx_baseline)) + sum(pxx_clean(idx_powerline)) + sum(pxx_clean(idx_emg));
signal_power = sum(pxx_orig(idx_signal));
SNR_before = 10 * log10(signal_power / noise_power_before);
SNR_after = 10 * log10(signal_power / noise_power_after);
SNR_improvement = SNR_after - SNR_before;
fprintf('\n========== SNR Report ==========\n');
fprintf('SNR Before Filtering : %.2f dB\n', SNR_before);
fprintf('SNR After  Filtering : %.2f dB\n', SNR_after);
fprintf('SNR Improvement      : %.2f dB\n', SNR_improvement);
fprintf('=================================\n');

 %poles and zeros 
figure;
subplot(1,3,1);
zplane(b_hp, a_hp);
title('HPF Pole-Zero');
subplot(1,3,2);
zplane(b_notch, a_notch);
title('Notch Pole-Zero');
subplot(1,3,3);
zplane(b_lp, a_lp);
title('LPF Pole-Zero');

% Impulse Responses
N_imp = 50;
N_imp_notch = 200;

impulse_s = [1; zeros(N_imp-1,1)];
impulse_s_notch = [1; zeros(N_imp_notch-1,1)];

figure;

subplot(1,3,1);
stem(filter(b_hp, a_hp, impulse_s));
title('HPF Impulse Response');
xlabel('n');
ylabel('h[n]');
grid on;

subplot(1,3,2);
stem(filter(b_notch, a_notch, impulse_s_notch));
title('Notch Impulse Response');
xlabel('n');
ylabel('h[n]');
grid on;

subplot(1,3,3);
stem(filter(b_lp, a_lp, impulse_s));
title('LPF Impulse Response');
xlabel('n');
ylabel('h[n]');
grid on;

% Step Responses
step_s       = ones(N_imp,1);
step_s_notch = ones(N_imp_notch, 1);
figure;
subplot(1,3,1);
stem(filter(b_hp, a_hp, step_s));
title('HPF Step Response');
xlabel('n');
ylabel('s[n]');
grid on;

subplot(1,3,2);
stem(filter(b_notch, a_notch, step_s_notch));
title('Notch Step Response');
xlabel('n');
ylabel('s[n]');
grid on;

subplot(1,3,3);
stem(filter(b_lp, a_lp, step_s));
title('LPF Step Response');
xlabel('n');
ylabel('s[n]');
grid on;

% Frequency Responses
[h_hp, f_axis] = freqz(b_hp, a_hp, 1024, Fs);
[h_notch, ~] = freqz(b_notch, a_notch, 1024, Fs);
[h_lp, ~] = freqz(b_lp, a_lp, 1024, Fs);
figure;
subplot(1,3,1);
plot(f_axis, 20*log10(abs(h_hp)));
title('HPF Magnitude (freq response)');
xlabel('Hz');
ylabel('dB');
grid on;
subplot(1,3,2);
plot(f_axis, 20*log10(abs(h_notch)));
title('Notch Magnitude (freq response)');
xlabel('Hz');
ylabel('dB');
grid on;
subplot(1,3,3);
plot(f_axis, 20*log10(abs(h_lp)));
title('LPF Magnitude (freq response)');
xlabel('Hz');
ylabel('dB');
grid on;
figure;
subplot(1,3,1);
plot(f_axis,angle(h_hp)*180/pi);
title('HPF Phase (freq response)');
xlabel('Hz');
ylabel('Degrees');
grid on;
subplot(1,3,2);
plot(f_axis, angle(h_notch)*180/pi);
title('Notch Phase (freq response)');
xlabel('Hz');
ylabel('Degrees');
grid on;
subplot(1,3,3);
plot(f_axis, angle(h_lp)*180/pi);
title('LPF Phase (freq response)');
xlabel('Hz');
ylabel('Degrees');
grid on;
close all;
clear;
clc;
data = readtable('C:/Users/Mega Store/Desktop/106_signals.csv');
ecg = data.MLII;
Fs = 360;
N = 3600;
ecg_segment = ecg(1:N);
t_segment = (0:N-1) / Fs;

Wp1 = (0.5/Fs) * 2 * pi;
Ws1 = (0.3/Fs) * 2 * pi;
Ap1 = 1;
As1 = 20;
OmegaP1 = 2 * Fs * tan(Wp1/2);
OmegaS1 = 2 * Fs * tan(Ws1/2);
[N1, ~] = buttord(OmegaP1, OmegaS1, Ap1, As1, 's');

Wp2 = ([49 51] / Fs) * 2 * pi;
Ws2 = ([49.5 50.5] / Fs) * 2 * pi;
Ap2 = 1;
As2 = 20;
OmegaP2 = 2 * Fs * tan(Wp2/2);
OmegaS2 = 2 * Fs * tan(Ws2/2);
[N2, ~] = buttord(OmegaP2, OmegaS2, Ap2, As2, 's');

Wp3 = (100 / Fs) * 2 * pi;
Ws3 = (120 / Fs) * 2 * pi;
Ap3 = 1;
As3 = 30;
OmegaP3 = 2 * Fs * tan(Wp3/2);
OmegaS3 = 2 * Fs * tan(Ws3/2);
[N3, ~] = buttord(OmegaP3, OmegaS3, Ap3, As3, 's');

[Z1, P1, K1] = buttap(N1);
[Z2, P2, K2] = buttap(N2);
[Z3, P3, K3] = buttap(N3);

[num1, den1] = zp2tf(Z1, P1, K1);
[num2, den2] = zp2tf(Z2, P2, K2);
[num3, den3] = zp2tf(Z3, P3, K3);

[num1, den1] = lp2hp(num1, den1, OmegaP1);

Wo2 = sqrt(OmegaP2(1) * OmegaP2(2));
Bw2 = OmegaP2(2) - OmegaP2(1);
[num2, den2] = lp2bs(num2, den2, Wo2, Bw2);

[num3, den3] = lp2lp(num3, den3, OmegaP3);

[bd1, ad1] = bilinear(num1, den1, Fs);
[bd2, ad2] = bilinear(num2, den2, Fs);
[bd3, ad3] = bilinear(num3, den3, Fs);

ecg_hp    = filtfilt(bd1, ad1, ecg_segment);
ecg_notch = filtfilt(bd2, ad2, ecg_hp);
ecg_clean = filtfilt(bd3, ad3, ecg_notch);

figure;
[pxx_orig,  f] = pwelch(ecg_segment, [], [], [], Fs);
[pxx_clean, ~] = pwelch(ecg_clean,   [], [], [], Fs);
plot(f, 10*log10(pxx_orig),  'r', 'LineWidth', 1.2); hold on;
plot(f, 10*log10(pxx_clean), 'b', 'LineWidth', 1.2);
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
colorbar;

subplot(2,1,2);
spectrogram(ecg_clean, 256, 200, 256, Fs, 'yaxis');
title('Filtered ECG Spectrogram');
colorbar;

figure;
subplot(2,1,1);
plot(t_segment, ecg_segment, 'r');
title('Original ECG');
ylabel('Amplitude (mV)');
grid on;

subplot(2,1,2);
plot(t_segment, ecg_clean, 'b');
title('Filtered ECG (Butterworth)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
grid on;

idx_baseline = f < 0.5;
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
fprintf('Filter Orders: HPF=%d, Notch=%d, LPF=%d\n', N1, N2, N3);


figure;
subplot(1,3,1);
zplane(bd1, ad1);
title('HPF Pole-Zero');
subplot(1,3,2);
zplane(bd2, ad2);
title('Notch Pole-Zero');
subplot(1,3,3);
zplane(bd3, ad3);
title('LPF Pole-Zero');


N_imp = 50;
N_imp_notch = 100;
impulse = [1; zeros(N_imp-1,       1)];
impulse_notch = [1; zeros(N_imp_notch-1, 1)];

figure;
tiledlayout(1, 3, 'Padding', 'normal', 'TileSpacing', 'normal');
nexttile;stem(filter(bd1, ad1, impulse)); 
title('HPF Impulse Response');
xlabel('n');
ylabel('h[n]');
grid on;
nexttile;
stem(filter(bd2, ad2, impulse_notch));
title('Notch Impulse Response');
xlabel('n');
ylabel('h[n]');
grid on;
nexttile;
stem(filter(bd3, ad3, impulse));
title('LPF Impulse Response');
xlabel('n');
ylabel('h[n]');
grid on;
step_sig   = ones(N_imp,1);
step_notch = ones(N_imp_notch, 1);
figure;
tiledlayout(1, 3, 'Padding', 'normal', 'TileSpacing', 'normal');
nexttile; stem(filter(bd1, ad1, step_sig));
title('HPF Step Response');
xlabel('n');
ylabel('s[n]');
grid on;
nexttile;
stem(filter(bd2, ad2, step_notch));
title('Notch Step Response');
xlabel('n');
ylabel('s[n]');
grid on;
nexttile;
stem(filter(bd3, ad3, step_sig));
title('LPF Step Response');
xlabel('n');
ylabel('s[n]');
grid on;
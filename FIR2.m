close all;
clear;
clc;
 
%% Signal Setup
data = readtable(['106' ...
    '_signals.csv']);
ecg  = data.MLII;
Fs   = 360;
N    = 3600;
ecg_segment = ecg(1:N);
t_segment   = (0:N-1) / Fs;
 
%% Filter Design
 
% FILTER 1: High-pass (baseline wander < 0.5 Hz)
Fst1 = 0.05;
Fp1  = 0.5;
Ap1  = 1;
As1  = 30;
Hf1 = fdesign.highpass('Fst,Fp,Ast,Ap', Fst1, Fp1, As1, Ap1, Fs);
Hd1 = design(Hf1, 'equiripple');
b1  = Hd1.Numerator;
a1  = 1;
 
% FILTER 2: FIR Notch (powerline at 50 Hz)
Fp1_n  = 47;
Fst1_n = 49;
Fst2_n = 51;
Fp2_n  = 53;
Ap2    = 0.5;
As2    = 40;
Hf2 = fdesign.bandstop('Fp1,Fst1,Fst2,Fp2,Ap1,Ast,Ap2', ...
      Fp1_n, Fst1_n, Fst2_n, Fp2_n, Ap2, As2, Ap2, Fs);
Hd2 = design(Hf2, 'equiripple');
b2  = Hd2.Numerator;
a2  = 1;
 
% FILTER 3: Low-pass (muscle noise > 100 Hz)
Fp3  = 100;
Fst3 = 110;
Ap3  = 0.5;
As3  = 40;
Hf3 = fdesign.lowpass('Fp,Fst,Ap,Ast', Fp3, Fst3, Ap3, As3, Fs);
Hd3 = design(Hf3, 'equiripple');
b3  = Hd3.Numerator;
a3  = 1;
 
%% Apply Filters
ecg_hp    = filter(b1, a1, double(ecg_segment));
ecg_notch = filter(b2, a2, ecg_hp);
ecg_clean = filter(b3, a3, ecg_notch);
 
%% Figure 1: PSD Comparison
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
 
%% Figure 2: Spectrograms
figure;
subplot(2,1,1);
spectrogram(ecg_segment, 256, 200, 256, Fs, 'yaxis');
title('Original ECG Spectrogram');
colorbar;
subplot(2,1,2);
spectrogram(ecg_clean, 256, 200, 256, Fs, 'yaxis');
title('Filtered ECG Spectrogram');
colorbar;
 
%% Figure 3: Time Domain
figure;
subplot(2,1,1);
plot(t_segment, ecg_segment, 'r');
title('Original ECG');
ylabel('Amplitude (mV)');
grid on;
subplot(2,1,2);
plot(t_segment, ecg_clean, 'b');
title('Filtered ECG (FIR Equiripple)');
xlabel('Time (s)');
ylabel('Amplitude (mV)');
grid on;
 
%% Figure 4: Pole-Zero Plots
figure;
subplot(1,3,1); zplane(b1, a1); title('HPF Pole-Zero');
subplot(1,3,2); zplane(b2, a2); title('Notch Pole-Zero');
subplot(1,3,3); zplane(b3, a3); title('LPF Pole-Zero');
 
%% Figure 5: Impulse Responses
N_imp       = 50;
N_imp_notch = 200;
impulse_s       = [1; zeros(N_imp-1, 1)];
impulse_s_notch = [1; zeros(N_imp_notch-1, 1)];
 
figure;
tiledlayout(1, 3, 'Padding', 'loose', 'TileSpacing', 'loose');
nexttile; stem(filter(b1, a1, impulse_s));
title('HPF Impulse Response'); xlabel('n'); ylabel('h[n]'); grid on;
nexttile; stem(filter(b2, a2, impulse_s_notch));
title('Notch Impulse Response'); xlabel('n'); ylabel('h[n]'); grid on;
nexttile; stem(filter(b3, a3, impulse_s));
title('LPF Impulse Response'); xlabel('n'); ylabel('h[n]'); grid on;
 
%% Figure 6: Step Responses
step_s       = ones(N_imp, 1);
step_s_notch = ones(N_imp_notch, 1);
 
figure;
tiledlayout(1, 3, 'Padding', 'loose', 'TileSpacing', 'loose');
nexttile; stem(filter(b1, a1, step_s));
title('HPF Step Response'); xlabel('n'); ylabel('s[n]'); grid on;
nexttile; stem(filter(b2, a2, step_s_notch));
title('Notch Step Response'); xlabel('n'); ylabel('s[n]'); grid on;
nexttile; stem(filter(b3, a3, step_s));
title('LPF Step Response'); xlabel('n'); ylabel('s[n]'); grid on;
 
%% Frequency Response Computation
[H1, w1] = freqz(b1, a1, 1024, Fs);
[H2, w2] = freqz(b2, a2, 1024, Fs);
[H3, w3] = freqz(b3, a3, 1024, Fs);
 
%% Figure 7: Magnitude Responses
figure;
tiledlayout(1, 3, 'Padding', 'loose', 'TileSpacing', 'loose');
nexttile; plot(w1, 20*log10(abs(H1)));
title('HPF Magnitude'); xlabel('Hz'); ylabel('dB'); grid on;
nexttile; plot(w2, 20*log10(abs(H2)));
title('Notch Magnitude'); xlabel('Hz'); ylabel('dB'); grid on;
nexttile; plot(w3, 20*log10(abs(H3)));
title('LPF Magnitude'); xlabel('Hz'); ylabel('dB'); grid on;
 
%% Figure 8: Phase Responses
figure;
tiledlayout(1, 3, 'Padding', 'loose', 'TileSpacing', 'loose');
nexttile; plot(w1, angle(H1)*(180/pi));
title('HPF Phase'); xlabel('Hz'); ylabel('Degrees'); grid on;
nexttile; plot(w2, angle(H2)*(180/pi));
title('Notch Phase'); xlabel('Hz'); ylabel('Degrees'); grid on;
nexttile; plot(w3, angle(H3)*(180/pi));
title('LPF Phase'); xlabel('Hz'); ylabel('Degrees'); grid on;
 
%% Figure 9: Group Delay
[gd1, wg1] = grpdelay(b1, a1, 1024, Fs);
[gd2, wg2] = grpdelay(b2, a2, 1024, Fs);
[gd3, wg3] = grpdelay(b3, a3, 1024, Fs);
 
figure;
plot(wg1, gd1/Fs * 1000); hold on;
plot(wg2, gd2/Fs * 1000);
plot(wg3, gd3/Fs * 1000);
xlabel('Frequency (Hz)');
ylabel('Group Delay (ms)');
title('Group Delay — FIR Filters');
legend('HPF', 'Notch', 'LPF');
xlim([0 150]);
grid on;
 
%% SNR Report
idx_baseline  = f < 0.5;
idx_powerline = f >= 49  & f <= 51;
idx_emg       = f >= 100 & f <= 150;
idx_signal    = f >= 0.5 & f <= 100;
 
noise_power_before = sum(pxx_orig(idx_baseline))  + ...
                     sum(pxx_orig(idx_powerline))  + ...
                     sum(pxx_orig(idx_emg));
noise_power_after  = sum(pxx_clean(idx_baseline))  + ...
                     sum(pxx_clean(idx_powerline))  + ...
                     sum(pxx_clean(idx_emg));
signal_power = sum(pxx_orig(idx_signal));
 
SNR_before      = 10 * log10(signal_power / noise_power_before);
SNR_after       = 10 * log10(signal_power / noise_power_after);
SNR_improvement = SNR_after - SNR_before;
 
fprintf('\n========== SNR Report ==========\n');
fprintf('SNR Before Filtering : %.2f dB\n', SNR_before);
fprintf('SNR After  Filtering : %.2f dB\n', SNR_after);
fprintf('SNR Improvement      : %.2f dB\n', SNR_improvement);
fprintf('=================================\n');
fprintf('Filter Orders: HPF=%d, Notch=%d, LPF=%d\n', ...
        order(Hd1), order(Hd2), order(Hd3));
 

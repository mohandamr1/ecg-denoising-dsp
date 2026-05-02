%% Signal Setup

data = readtable('100_signals.csv');

ecg = data.MLII;
Fs  = 360;  
N = 3600;


ecg_segment = ecg(1:N);  
t_segment = (0:N-1) / Fs;
frq = (-N/2 : N/2 - 1) * Fs/N;
%% Chebyshev type II design (0.5-100)Hz

% FILTER 1: Baseline HPF < 0.5

Wp1 = (0.5/Fs) * 2 * pi;   % passband edge (normalized)
Ws1 = (0.3/Fs) * 2 * pi;   % stopband edge
Rp1 = 1;           % passband ripple (dB)
Rs1 = 30;          % stopband attenuation (dB)

OmegaP1 = 2 * Fs * tan(Wp1/2);
OmegaS1 = 2 * Fs * tan(Ws1/2);

[N1, Wn1] = cheb2ord(OmegaP1, OmegaS1, Rp1, Rs1,'s');

% FILTER 2: Notch at 50 Hz
Wp2 = ([49 51] / Fs)* 2 * pi;   % passband (just outside 50 Hz)
Ws2 = ([49.5 50.5] / Fs) * 2 * pi; % stopband (tight around 50 Hz)
Rp2 = 1;
Rs2 = 20;


OmegaP2 = 2 * Fs * tan(Wp2/2);
OmegaS2 = 2 * Fs * tan(Ws2/2);

[N2, Wn2] = cheb2ord(OmegaP2, OmegaS2, Rp2, Rs2,'s');

% FILTER 3: Low-pass (muscle noise) (20-150Hz)
Wp3 = (100 / Fs) * 2 * pi;   % passband edge
Ws3 = (120 / Fs) * 2 * pi;   % stopband edge
Rp3 = 1;
Rs3 = 40;

OmegaP3 = 2 * Fs * tan(Wp3/2);
OmegaS3 = 2 * Fs * tan(Ws3/2);

[N3, Wn3] = cheb2ord(OmegaP3, OmegaS3, Rp3, Rs3,'s');

[Z1,P1,K1] = cheb2ap(N1,Rp1);
[Z2,P2,K2] = cheb2ap(N2,Rp2);
[Z3,P3,K3] = cheb2ap(N3,Rp3);


[num1,den1] = zp2tf(Z1,P1,K1);
[num2,den2] = zp2tf(Z2,P2,K2);
[num3,den3] = zp2tf(Z3,P3,K3);

[num1,den1] = lp2hp(num1, den1, OmegaP1);   % LP to HP
[bd1,ad1]   = bilinear(num1, den1, Fs); 


Wo2 = sqrt(OmegaP2(1) * OmegaP2(2));        % center frequency
Bw2 = OmegaP2(2) - OmegaP2(1);             % bandwidth

[num2,den2] = lp2bs(num2, den2, Wo2, Bw2);  % LP to BS
[bd2,ad2]   = bilinear(num2, den2, Fs); 


[num3,den3] = lp2lp(num3, den3, OmegaP3);  
[bd3,ad3]   = bilinear(num3, den3, Fs);  

ecg_hp    = filtfilt(bd1, ad1, ecg_segment);
ecg_notch = filtfilt(bd2, ad2, ecg_hp);
ecg_clean = filtfilt(bd3, ad3, ecg_notch);
disp(N1);
disp(N2);
disp(N3);
figure;
[pxx_orig,  f] = pwelch(ecg_segment, [], [], [], Fs);  % ← f not f_orig
[pxx_clean, ~] = pwelch(ecg_clean,   [], [], [], Fs);

plot(f, 10*log10(pxx_orig),  'r'); hold on;
plot(f, 10*log10(pxx_clean), 'b');
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title('PSD: Original vs Filtered ECG');
legend('Original', 'Filtered');
xlim([0 180]);
grid on;

% Spectrograms
figure;
subplot(2,1,1);
spectrogram(ecg_segment, 256, 200, 256, Fs, 'yaxis');
title('Original ECG Spectrogram');
subplot(2,1,2);
spectrogram(ecg_clean, 256, 200, 256, Fs, 'yaxis');
title('Filtered ECG Spectrogram');


figure;
subplot(2,1,1); plot(t_segment, ecg_segment); title('Original ECG');
subplot(2,1,2); plot(t_segment, ecg_clean);   title('Filtered ECG');
xlabel('Time (s)');


idx_baseline  = f < 0.5;
idx_powerline = f >= 49 & f <= 51;
idx_emg       = f >= 100 & f <= 150;
idx_signal    = f >= 0.5 & f <= 100;

noise_power_before = sum(pxx_orig(idx_baseline)) + ...
                     sum(pxx_orig(idx_powerline)) + ...
                     sum(pxx_orig(idx_emg));

noise_power_after  = sum(pxx_clean(idx_baseline)) + ...
                     sum(pxx_clean(idx_powerline)) + ...
                     sum(pxx_clean(idx_emg));

signal_power = sum(pxx_orig(idx_signal));

SNR_before = 10*log10(signal_power / noise_power_before);
SNR_after  = 10*log10(signal_power / noise_power_after);

fprintf('SNR Before:      %.2f dB\n', SNR_before);
fprintf('SNR After:       %.2f dB\n', SNR_after);
fprintf('SNR Improvement: %.2f dB\n', SNR_after - SNR_before);



figure;
subplot(1,3,1); zplane(bd1, ad1); title('HPF Pole-Zero');
subplot(1,3,2); zplane(bd2, ad2); title('Notch Pole-Zero');
subplot(1,3,3); zplane(bd3, ad3); title('LPF Pole-Zero');



N = 50;  
N_notch = 100;  
impulse = [1; zeros(N-1, 1)];
impulse_notch = [1; zeros(N_notch-1, 1)];


step_sig      = ones(N, 1);
step_notch    = ones(N_notch, 1);

figure;
tiledlayout(1,3, 'Padding', 'loose', 'TileSpacing', 'loose');

nexttile; stem(filter(bd1,ad1,impulse));       title('HPF Impulse');
nexttile; stem(filter(bd2,ad2,impulse_notch)); title('Notch Impulse');
nexttile; stem(filter(bd3,ad3,impulse));       title('LPF Impulse');

% Step Response
figure;
tiledlayout(1,3, 'Padding', 'loose', 'TileSpacing', 'loose');

nexttile; stem(filter(bd1,ad1,step_sig));      title('HPF Step');
nexttile; stem(filter(bd2,ad2,step_notch));    title('Notch Step');
nexttile; stem(filter(bd3,ad3,step_sig));      title('LPF Step');     title('LPF Step');



[H1, w1] = freqz(bd1, ad1, 1024, Fs);
[H2, w2] = freqz(bd2, ad2, 1024, Fs);
[H3, w3] = freqz(bd3, ad3, 1024, Fs);

% Magnitude Response
figure;
tiledlayout(1,3, 'Padding', 'loose', 'TileSpacing', 'loose');

nexttile; plot(w1, 20*log10(abs(H1))); 
title('HPF Magnitude'); xlabel('Hz'); ylabel('dB'); grid on;

nexttile; plot(w2, 20*log10(abs(H2))); 
title('Notch Magnitude'); xlabel('Hz'); ylabel('dB'); grid on;

nexttile; plot(w3, 20*log10(abs(H3))); 
title('LPF Magnitude'); xlabel('Hz'); ylabel('dB'); grid on;

% Phase Response
figure;
tiledlayout(1,3, 'Padding', 'loose', 'TileSpacing', 'loose');

nexttile; plot(w1, angle(H1)*(180/pi)); 
title('HPF Phase'); xlabel('Hz'); ylabel('Degrees'); grid on;

nexttile; plot(w2, angle(H2)*(180/pi)); 
title('Notch Phase'); xlabel('Hz'); ylabel('Degrees'); grid on;

nexttile; plot(w3, angle(H3)*(180/pi)); 
title('LPF Phase'); xlabel('Hz'); ylabel('Degrees'); grid on;
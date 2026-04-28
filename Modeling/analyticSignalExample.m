%% Analytic signal + negative frequencies + complex baseband demo
% Requires Signal Processing Toolbox for hilbert().
% If you don't have it, tell me and I'll give an FFT-only Hilbert version.

clear; close all; clc;

%% 1) Sampling and time axis
fs = 2000;                 % sample rate [Hz]
T  = 2.0;                  % duration [s]
t  = (0:1/fs:T-1/fs).';    % column vector
N  = numel(t);

%% 2) Build a "passband" signal x(t) = A(t) cos(2pi f_c t + phi(t))
fc = 250;                  % carrier [Hz]

% Slow amplitude envelope A(t): a smooth pulse (Gaussian)
A = exp(-((t - T/2)/(0.25)).^2);

% Slow phase term phi(t): make it a chirp-like phase (quadratic phase)
% A linear instantaneous frequency offset is produced by quadratic phase.
% phi(t) = 2*pi*(f_offset * t + 0.5*k*t^2) (baseband phase)
f0 = -40;                  % baseband start offset [Hz]
f1 =  40;                  % baseband end offset [Hz]
k  = (f1 - f0)/T;          % chirp rate [Hz/s]
phi = 2*pi*( f0*t + 0.5*k*t.^2 );

% Real passband signal
x = A .* cos(2*pi*fc*t + phi);

%% 3) Analytic signal via Hilbert transform: x_a(t) = x(t) + i*Hilbert{x(t)}
xa = hilbert(x);  % complex analytic signal

% Envelope and phase of analytic signal (passband envelope, passband phase)
env_xa   = abs(xa);
phase_xa = unwrap(angle(xa));

%% 4) Show "negative frequencies" vs "one-sided" analytic spectrum
% FFT helpers
f = (-N/2:N/2-1).' * (fs/N);

X  = fftshift(fft(x));
XA = fftshift(fft(xa));

% Magnitude spectra (in dB, normalized)
Xmag  = 20*log10( abs(X)  / max(abs(X))  + 1e-12 );
XAmag = 20*log10( abs(XA) / max(abs(XA)) + 1e-12 );

figure('Name','Negative frequencies & analytic spectrum');
subplot(2,1,1);
plot(f, Xmag, 'LineWidth', 1.2); grid on;
xlabel('Frequency (Hz)'); ylabel('|X(f)| (dB)');
title('Real signal spectrum: symmetric about 0 (positive + negative frequencies)');
xlim([-600 600]);

subplot(2,1,2);
plot(f, XAmag, 'LineWidth', 1.2); grid on;
xlabel('Frequency (Hz)'); ylabel('|X_a(f)| (dB)');
title('Analytic signal spectrum: negative frequencies suppressed');
xlim([-600 600]);

%% 5) Mix down to complex baseband: s(t) = x_a(t) * exp(-i2pi f_c t)
s_raw = xa .* exp(-1i*2*pi*fc*t);

% In ideal math, this already places content near 0 Hz, but you still
% typically low-pass filter to remove any residual high-frequency terms.
% We'll apply a simple lowpass to mimic receiver filtering.
s = lowpass(s_raw, 120, fs);   % cutoff > max baseband sweep (~40 Hz)

env_s   = abs(s);
phase_s = unwrap(angle(s));

%% 6) Compare envelope and phase to what we injected (A(t), phi(t))
% The baseband phase should match phi(t) up to an arbitrary constant offset.
% We can align by subtracting the initial phase difference.
phase_s_aligned = phase_s - (phase_s(1) - phi(1));

figure('Name','Complex baseband: envelope and phase');
subplot(3,1,1);
plot(t, x, 'LineWidth', 1.0); grid on;
xlabel('Time (s)'); ylabel('x(t)');
title('Real passband signal x(t)');

subplot(3,1,2);
plot(t, env_s, 'LineWidth', 1.2); hold on;
plot(t, A, '--', 'LineWidth', 1.2);
grid on; xlabel('Time (s)'); ylabel('Amplitude');
title('Baseband envelope |s(t)| tracks A(t)');
legend('|s(t)|','A(t)');

subplot(3,1,3);
plot(t, phase_s_aligned, 'LineWidth', 1.2); hold on;
plot(t, phi, '--', 'LineWidth', 1.2);
grid on; xlabel('Time (s)'); ylabel('Phase (rad)');
title('Baseband phase angle(s(t)) tracks \phi(t) (up to constant offset)');
legend('unwrap(angle(s))','\phi(t)');

%% 7) Show baseband spectrum centered near 0 Hz
S = fftshift(fft(s));
Smag = 20*log10( abs(S)/max(abs(S)) + 1e-12 );

figure('Name','Baseband spectrum near DC');
plot(f, Smag, 'LineWidth', 1.2); grid on;
xlabel('Frequency (Hz)'); ylabel('|S(f)| (dB)');
title('Complex baseband spectrum: content shifted near 0 Hz (DC)');
xlim([-200 200]);

%% 8) Optional: show I/Q trajectories (phasor view)
figure('Name','I/Q phasor view');
plot(real(s), imag(s), 'LineWidth', 1.0); grid on; axis equal;
xlabel('I = Re\{s(t)\}'); ylabel('Q = Im\{s(t)\}');
title('Complex baseband in I/Q plane (phasor trajectory)');

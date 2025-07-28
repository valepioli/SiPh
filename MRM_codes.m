%% ========================
%  MRM code 1
% =========================
clear all; % Clear all variables from the workspace

group_index = 4.07849;
% Group index of the waveguide mode (dimensionless)

group_index_vvar = group_index * (2.0e-5);
% Group index with a small variation (e.g., due to tuning or process variation)

Radius = 86.6e-6;
% Ring radius in meters
L = 2 * pi * Radius;
% Ring Circumference (m)

c = 3e8;
% Speed of Light in vacuum (m/s), positive value

te = 2.0465 * 1e-11;
% External coupling time constant (s)

ti = 1.282 * 1e-11;
% Internal loss time constant (s)

n = 855;
% Mode number of the resonator

w0 = 2 * pi * c * n / (group_index * L);
% Resonance angular frequency (rad/s)

w0_lambda = 2 * pi * c / w0;
% Resonance wavelength corresponding to w0 (m)

w0_var = 2 * pi * c * n / (group_index_vvar * L);
% Resonance angular frequency with varied group index (rad/s)

w0_lambda_vvar = 2 * pi * c / w0_var;
% Resonance wavelength corresponding to w0_var (m)

w0_shift = w0_lambda_vvar - w0_lambda;
% Wavelength shift due to group index variation (m)

Dropshift = w0_shift * 1e12;
% Convert wavelength shift to picometers (pm)

w0_m1 = 2 * pi * c * (n + 1) / (group_index * L);
% Angular frequency of the next mode (n + 1)

w0_lambda_m1 = 2 * pi * c / w0_m1;
% Wavelength corresponding to the next mode (m)

FSR = w0_lambda_m1 - w0_lambda;
% Free spectral range (FSR) in the wavelength domain (m)


%% ========================
%  MRM code 2
% =========================
clear all; clc;

% == Parameter Settings ==
w0_lambda_target = 1559.5e-9;     % Target resonance wavelength at -IV (in meters)
group_index_0 = 4.07149;          % Group index at 0V
group_index_4V = 4.07849;         % Group index at -4V
group_index_diff = 6.2e-3;        % Δn_g = n_g(0V) - n_g(4V)
group_index_vvar = group_index_0 - group_index_diff; % Group index at -IV
Radius = 86.6e-6;                 % Ring radius (in meters)
L = 2 * pi * Radius;              % Ring circumference (in meters)

% == Find the closest mode number m such that w_lambda is close to target ==
m_try = round(group_index_vvar * L / w0_lambda_target);
lambda_test = group_index_vvar * L / m_try;

while lambda_test > w0_lambda_target
    m_try = m_try + 1;
    lambda_test = group_index_vvar * L / m_try;
end

lambda_down = group_index_vvar * L / (m_try - 1);
if abs(lambda_down - w0_lambda_target) < abs(lambda_test - w0_lambda_target)
    m_try = m_try - 1;
    lambda_vvar = lambda_down;
else
    lambda_vvar = lambda_test;
end

lambda_0V = group_index_0 * L / m_try;  
Dropshift = (lambda_vvar - lambda_0V) * 1e12;

% == Print Results ==
fprintf('Final mode number = %d\n', m_try);
fprintf('Wavelength at 0V = %.4f nm\n', lambda_0V * 1e9);
fprintf('Wavelength at -IV = %.4f nm\n', lambda_vvar * 1e9);
fprintf('Wavelength shift = %.4f pm\n', Dropshift);

% == Simulate Spectrum ==
lambda = linspace(lambda_0V - 2.0e-9, lambda_0V + 2.0e-9, 2000); % Wavelength sweep range (12nm)
lambda_0 = 2e-12;                                                % Full-width at half maximum ~80pm
delta_lambda = lambda - lambda_0V;
transmission = 1 - exp(-(delta_lambda / (lambda_0V/2)).^2);      % Gaussian shape transmission

P_in = 1e-3;                        % Input power (1 mW)
P_out = P_in * transmission;       % Output power (in Watts)
P_dB = 10 * log10(P_out / P_in);   % Relative power in dB

% == Wavelengths corresponding to -3 dB, -5 dB, -6 dB levels ==
target_dB = [-3, -5, -6];
thresholds = 10.^(target_dB / 10); % Convert dB to power ratio
lambda_marks = zeros(length(target_dB), 2); % Store [left, right] wavelengths

for k = 1:length(target_dB)
    idx_left = find(P_out < thresholds(k) * P_in, 1, 'first');
    idx_right = idx_left + find(P_out(idx_left:end) > thresholds(k) * P_in, 1, 'first') - 1;
    lambda_marks(k, :) = [lambda(idx_left), lambda(idx_right)];
end

% == Plotting ==
figure;

% --- Plot in Watt scale ---
subplot(2,1,1);
plot(lambda * 1e9, P_out, 'b', 'LineWidth', 1.5); hold on;
xlabel('Wavelength (nm)');
ylabel('Power (W)');
title('Spectrum (Watt scale)');
grid on;
ylim([0, 1e-3]);
yticks(0:0.2e-3:1e-3);

% --- Plot in dB scale ---
subplot(2,1,2);
plot(lambda * 1e9, P_dB, 'r', 'LineWidth', 1.5); hold on;
xlabel('Wavelength (nm)');
ylabel('Transmission (dB)');
title('Spectrum (dB scale)');
grid on;
ylim([-30, 0]);
yticks(-30:20:0);

% == Mark -3/-5/-6 dB crossover wavelengths ==
color_map = [0 0 0; 0 0.5 0; 0.5 0 0.5]; % Black, Green, Purple
yl = ylim;

for k = 1:length(target_dB)
    xL = lambda_marks(k, 1) * 1e9;
    xR = lambda_marks(k, 2) * 1e9;
    yval = target_dB(k);

    % Vertical lines
    line([xL xL], yl, 'color', color_map(k,:), 'LineStyle', '--', 'LineWidth', 1.2);
    line([xR xR], yl, 'color', color_map(k,:), 'LineStyle', '--', 'LineWidth', 1.2);

    % Text annotations
    text(xL, yval, sprintf('%.4f nm', xL), 'color', color_map(k,:), ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    text(xR, yval, sprintf('%.4f nm', xR), 'color', color_map(k,:), ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end

legend('spectrum', '-3dB', '-5dB', '-6dB', 'Location', 'southwest');

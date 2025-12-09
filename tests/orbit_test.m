% -------------------------------------------------------------------------
% SIMULINK SIMULATION RUNNER: ISS ORBIT MAG FIELD VERIFICATION
% -------------------------------------------------------------------------
clearvars; clc; close all;

% Loading paths
addpath('../src');
addpath('../data');
addpath('../plots'); 

% Loading coefficients
load('../data/igrf14coeffs.mat', 'C_nm');

% Simulation Parameters

% Physical Constants
Re = 6378137;           % Earth Radius (m)
mu = 3.986004418e14;    % Std Gravitational Parameter

% ISS Orbit 
h_p = 415000;           % Altitude ~415 km
inc = 51.6;             % Inclination (deg)
e   = 0.0006;           % Eccentricity (Near circular)
omega = 0;              % Argument of Periapsis
RAAN  = 0;              % Right Ascension
nu    = 0;              % True Anomaly start

% Derived Elements
a = Re + h_p;           % Semi-major axis
T_period = 2*pi * sqrt(a^3 / mu);

% Simulation Time
num_orbits = 2;
t_stop = num_orbits*T_period; % 2 full orbits

% Epoch & Model Degree
dec_year = 2025.0;      % Decimal Year 
degree   = 13;          % Harmonic Degree for Custom Model

model_name = 'orbit_simulation';
try
    out = sim(model_name);
catch ME
    error('Failed to run Simulink model: %s. \nMake sure "%s.slx" is in the path or current folder.', ME.message, model_name);
end

% Data extraction
extract_data = @(ts) squeeze(ts.Data); 

t = out.B_custom.Time;

% B vectors 
B_cust = extract_data(out.B_custom);      
B_ref  = extract_data(out.B_igrfmagm);    

% Magnitudes
F_cust = extract_data(out.F_custom);
F_ref  = extract_data(out.F_igrfmagm);

% Dimentions Handling
if size(B_cust, 2) ~= 3 && size(B_cust, 1) == 3
    B_cust = B_cust';
end
if size(B_ref, 2) ~= 3 && size(B_ref, 1) == 3
    B_ref = B_ref';
end

n_min = min(size(B_cust, 1), size(B_ref, 1));
B_cust = B_cust(1:n_min, :);
B_ref  = B_ref(1:n_min, :);
t = t(1:n_min);

% Error Computation
err_vec = B_cust - B_ref;

% Magnitude of the Error Vector
err_mag = vecnorm(err_vec, 2, 2); 

% Plots
fig1 = figure('Color', 'w', 'Position', [100 100 1000 800]);
sgtitle(['IGRF-14 Verification: Custom vs Toolbox (N=' num2str(degree) ') (ISS Orbit)']);

t_min = t/60;

subplot(3,1,1);
plot(t_min, B_ref(:,1), 'k', 'LineWidth', 2); hold on;
plot(t_min, B_cust(:,1), 'r--', 'LineWidth', 1.5);
ylabel('North (Bx) [nT]'); grid on; legend('igrfmagm', 'Custom Model');
title('Magnetic Vector Components');

subplot(3,1,2);
plot(t_min, B_ref(:,2), 'k', 'LineWidth', 2); hold on;
plot(t_min, B_cust(:,2), 'r--', 'LineWidth', 1.5);
ylabel('East (By) [nT]'); grid on;

subplot(3,1,3);
plot(t_min, B_ref(:,3), 'k', 'LineWidth', 2); hold on;
plot(t_min, B_cust(:,3), 'r--', 'LineWidth', 1.5);
ylabel('Down (Bz) [nT]'); xlabel('Time [min]'); grid on;

saveas(fig1, '../plots/simulink_vector_comparison.png');

% Error Analysis
fig2 = figure('Color', 'w', 'Position', [150 150 800 400]);
plot(t_min, err_mag, 'b', 'LineWidth', 1.5);
grid on;
xlabel('Time [min]'); ylabel('Error Magnitude [nT]');
title('Absolute Vector Error (|B_{custom} - B_{ref}|)');

% Statistics
mean_err = mean(err_mag);
max_err  = max(err_mag);
yline(mean_err, '--r', ['Mean: ' num2str(mean_err, '%.3f') ' nT']);
subtitle(['Max Error: ' num2str(max_err, '%.3f') ' nT']);

saveas(fig2, '../plots/simulink_error_analysis.png');

fprintf('Plots saved to ../plots/\n');
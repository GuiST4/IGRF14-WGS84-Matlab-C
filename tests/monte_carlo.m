% Monte Carlo Accuracy Test - Custom Model VS IGRF-13
clearvars; clc; close all;

% Loading paths
addpath('../src');
addpath('../data');
addpath('../tools');

% Load coefficients
load('../data/igrf14coeffs.mat', 'C_nm');

% Simulation Parameters
n_sim = 200;           % Number of random points to test
max_degree = 13;       % Maximum degree to test
year = 2025.0;         % Test epoch

% Random Geodetic Coordinates (LEO Orbit: 300-800km)
% Lat: -90 to 90, Lon: -180 to 180
rng(42); 
lats = -90 + 180 * rand(n_sim, 1);
lons = -180 + 360 * rand(n_sim, 1);
alts = (300 + 500 * rand(n_sim, 1))*1000;


rmse_results = zeros(max_degree, 1);

fprintf('Monte Carlo Simulation (n=%d)\n', n_sim);

% Aerospace toolbox is needed!!
if isempty(which('igrfmagm'))
    error('Aerospace Toolbox not found.');
end

% Loop from degree n = 0 to n = 13
for n = 1:max_degree

    errors_sq = 0;
    
    for i = 1:n_sim

        [~, F_custom] = magnetic_field(alts(i), lats(i), lons(i), year, C_nm, n);
        [~, ~, ~, ~, F_ref] = igrfmagm(alts(i), lats(i), lons(i), year);
        
        % Accumulate Squared Error
        errors_sq = errors_sq + (F_custom - F_ref)^2;
    end
    
    % Compute RMSE (Root Mean Square Error)
    rmse_results(n) = sqrt(errors_sq / n_sim);
    fprintf('Degree N=%2d | RMSE: %.4f nT\n', n, rmse_results(n));

end

% Plot Results
figure('Color', 'w', 'Name', 'RMSE Analysis');
semilogy(1:max_degree, rmse_results, '-o', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on;
title(['RMSE | Custom Model VS igrfmagm | (Monte Carlo: n = ' num2str(n_sim) ')']);
xlabel('Spherical Harmonic Degree (N)');
ylabel('RMSE Magnitude (nT) [Log Scale]');
xlim([1 max_degree]);
xticks(1:max_degree);

% Sensor noise (approx 200nT for cheap sensors)
yline(200, '--b', 'Sensor Noise (~200nT)');
legend('Model RMSE', 'Sensor Noise','Location','southeast');

plot_dir = '../plots';
output_filename = fullfile(plot_dir, 'rmse_convergence.png');
saveas(gcf, output_filename);
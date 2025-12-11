# IGRF-14 Geomagnetic Field Model (WGS84) - MATLAB & C

A high-fidelity implementation of the International Geomagnetic Reference Field (IGRF-14) model for Low Earth Orbit (LEO) applications. 

This core module is designed as a reference verification tool (**MATLAB**) and a flight-ready driver (**Embedded C**) for CubeSat Attitude Determination and Control Systems (ADCS).

## Features

- **High Fidelity:** Implements IGRF-14 generation up to degree/order 13.
- **Geodetic Accuracy:** Includes a rigorous WGS84 Geodetic-to-Geocentric coordinate transformation (accounting for Earth's oblateness).
- **Embedded C Driver:** Static memory allocation (no `malloc`), single-precision arithmetic (`float`), and zero external dependencies.
- **Secular Variation:** Linear time interpolation for precise field estimation between epoch years.

## Usage

### 1. MATLAB (Analysis & Verification)

The MATLAB implementation is located in [src/](src/) and can be run as follows:

```matlab
% Load coefficients
load('../data/igrf14coeffs.mat', 'C_nm');

% Define position 
height = 400000;    % 400 km
lat    = 38.7223;   % Lisbon Latitude
lon    = -9.1393;   % Lisbon Longitude
year   = 2025.5;    % Date
N      = 13;        % Expansion order

% Compute Field
[B_ned, F] = magnetic_field(height, lat, lon, year, C_nm, N);

disp(['Total Intensity: ', num2str(F), ' nT']);
```

### Embedded C (Flight Software)

The C library is located in [src_c/](src_c/). It uses auto-generated static arrays for coefficients.
```c
#include "src_c/igrf.h"

// Define Inputs
igrf_input_t in = {
    .lat_deg = 38.7223f,
    .lon_deg = -9.1393f,
    .alt_m   = 400000.0f,
    .year    = 2025.5f
};

igrf_output_t out;

// Run Model
igrf_compute(&in, &out);

// Result is in out.B_north, out.B_east, out.B_down
```

## Verification Results

### 1. Vector Comparison (Custom vs MATLAB Aerospace Toolbox)
The following plot compares the North, East, and Down magnetic vector components over 2 orbits.
<p align="center">
  <img src="plots/simulink_vector_comparison.png" width="80%">
</p>

### 2. RMSE Convergence
A Monte Carlo analysis showing the Root Mean Square Error (RMSE) of the custom model converging to the MATLAB reference as the spherical harmonic degree (N) increases.
<p align="center">
  <img src="plots/rmse_convergence.png" width="80%">
</p>

## Disclaimer
This model is validated against the Aerospace Toolbox. However, if a sign error in this code causes the satellite to spin up like a Beyblade instead of stabilizing:
1. I am not liable.
2. You win a high-speed centrifuge record.

## License
MIT License

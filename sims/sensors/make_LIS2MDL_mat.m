% make_LIS2MDL_mat  Creates LIS2MDL.mat (magnetometer chip parameters).
% Numbers come from the LIS2MDL datasheet (DS12144 Rev 6). MATLAB wants microtesla (uT):
%   1 gauss = 100 uT,   1 mG (milligauss) = 0.1 uT
% Earth's field is NOT stored here; sim_driver computes it with wrldmagm.

mag = struct();   % start with an empty container so an old 'mag' can't interfere

% Biggest field it can read: +/-49.152 gauss (Table 2, FS)
mag.MeasurementRange = 4915.2;                 % uT

% Smallest step it can read: 1.5 mG per count (Table 2, So)
mag.Resolution = 0.15;                         % uT per count

% Noise: 4.5 mG RMS = 0.45 uT (Table 9, high-res, no LPF/offset cancellation).
% WARNING: assumes 100 Hz. MATLAB wants density, not total noise:
% total = density * sqrt(rate/2). If the sample rate changes, change the 100 below.
mag.NoiseDensity = [1 1 1] * 0.45 / sqrt(100/2);    % = 0.0636 uT/sqrt(Hz)

% Sensitivity drift with temperature: +/-0.03 %/degC (Table 2, TCSo)
mag.TemperatureScaleFactor = [1 1 1] * 0.03;        % %/degC

% Offset drift with temperature: up to 0.3 mG per degC (Table 2, TCOff)
mag.TemperatureBias = [1 1 1] * 0.03;               % uT per degC

save("LIS2MDL.mat", "mag");
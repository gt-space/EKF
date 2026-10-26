classdef LIS2MDL < handle
    % LIS2MDL  Magnetometer sensor model (same pattern as ASM330LHBG1.m).
    %
    %   What it does: takes the rocket's TRUE attitude and Earth's TRUE
    %   magnetic field, and returns what our magnetometer would READ:
    %   the field seen from the rocket's body frame, plus realistic noise.
    %
    %   Usage:
    %       mag = LIS2MDL(100, "LIS2MDL.mat");             % 100 Hz = RECO's setting
    %       reading = mag.measure(simData.quat(:,k), B);    % 3x1, microtesla
    %   where B = Earth's field at that instant, NED, uT (from sim_driver)

    properties (SetAccess = private)
        name = "LIS2MDL"
        sample_rate      % Hz
        num_measurements % number of measurements
        sensor_category
        file_name = ""
        imu
    end

    methods
        function obj = LIS2MDL(sample_rate, arg2)

            obj.sample_rate = sample_rate;
            obj.num_measurements = 3;   % Bx, By, Bz

            % MATLAB has no magnetometer-only imuSensor; 'accel-mag' is the
            % smallest option with a magnetometer. We ignore the accel output.
            obj.imu = imuSensor( ...
                'accel-mag', ...
                'SampleRate', sample_rate);

            if nargin ~= 2
                error("LIS2MDL:InvalidArguments", ...
                    "Usage: LIS2MDL(sampleRate, fileName) or LIS2MDL(sampleRate, magParams).");
            end

            if isstruct(arg2)
                % Usage: LIS2MDL(sampleRate, magParams)
                magParams = arg2;
            else
                % Usage: LIS2MDL(sampleRate, fileName)
                obj.file_name = string(arg2);
                loaded = load(obj.file_name);

                if ~isfield(loaded, 'mag')
                    error("LIS2MDL:InvalidFile", ...
                        "File must contain a 'mag' struct.");
                end

                magParams = loaded.mag;
            end

            % Build a fresh parameter object and copy each datasheet number in.
            % (Assigning the whole object at once is what makes the settings stick.)
            mp = magparams;
            fields = fieldnames(magParams);
            for i = 1:numel(fields)
                if ~isprop(mp, fields{i})
                    error("LIS2MDL:InvalidParameter", ...
                        "Unsupported sensor parameter: %s", fields{i});
                end
                mp.(fields{i}) = magParams.(fields{i});
            end
            obj.imu.Magnetometer = mp;

            obj.sensor_category = measurement_type.MAG;
        end

        function measurement = measure(obj, trueQuat, trueField)
            % trueQuat:  4x1 [w x y z], body -> NED (same as simData.quat)
            % trueField: 3x1 Earth's field at this instant, NED, uT
            %
            % imuSensor rotates the NED field into the body frame using the
            % orientation, then adds the noise/quantization from the .mat.
            obj.imu.MagneticField = trueField(:).';
            orientation = quaternion(trueQuat(:).');
            [~, measuredMag] = obj.imu(zeros(1,3), zeros(1,3), orientation);

            measurement = measuredMag.';   % 3x1, microtesla
        end
    end
end
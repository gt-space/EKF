classdef ASM330LHBG1 < handle

    properties (SetAccess = private)
        name = "ASM330LHBG1"
        sample_rate      % Hz
        num_measurements % number of measurements
        sensor_catgeory
        file_name = ""
        imu
    end

    methods
        function obj = ASM330LHBG1(sample_rate, arg2, arg3)

            obj.sample_rate = sample_rate;
            obj.num_measurements = 6;

            % Initialize the MATLAB IMU model.
            obj.imu = imuSensor( ...
                'accel-gyro', ...
                'SampleRate', sample_rate);

            if nargin == 3
                % Usage: ASM330LHBG1(sampleRate, accelParams, gyroParams)
                accelParams = arg2;
                gyroParams = arg3;

            elseif nargin == 2
                % Usage: ASM330LHBG1(sampleRate, fileName)
                obj.file_name = string(arg2);

                loaded = load(obj.file_name);

                if ~isfield(loaded, 'accel') || ...
                   ~isfield(loaded, 'gyro')
                    error("ASM330LHBG1:InvalidFile", ...
                        "File must contain 'accel' and 'gyro' structs.");
                end

                accelParams = loaded.accel;
                gyroParams = loaded.gyro;

            else
                error("ASM330LHBG1:InvalidArguments", ...
                    "Provide sampleRate and either a filename " + ...
                    "or accelerometer and gyroscope parameters.");
            end

            obj.configureSensor(obj.imu.Accelerometer, accelParams);
            obj.configureSensor(obj.imu.Gyroscope, gyroParams);
            obj.sensor_catgeory = measurement_type.IMU;
        end

        function measurement = measure(obj, trueAcceleration, trueAngularVelocity)

            [measuredAccel, measuredAngularVelocity] = obj.imu( ...
                trueAcceleration, trueAngularVelocity);

           measurement = [measuredAccel; measuredAngularVelocity];
        end
    end

    methods (Access = private)
        function configureSensor(~, sensor, sensorParams)

            fields = fieldnames(sensorParams);

            for i = 1:numel(fields)
                field = fields{i};

                if ~isprop(sensor, field)
                    error("ASM330LHBG1:InvalidParameter", ...
                        "Unsupported sensor parameter: %s", field);
                end

                sensor.(field) = sensorParams.(field);
            end
        end
    end
end
classdef sim_driver

    properties
        timestep
        sensors
        sensor_gather_intervals
        sensor_measurement_types
        main_filter
        backup_filter
        parachute_deployment_logic
        truth
        truth_hash
        truth_file_path
        output_dir_path
    end

    methods
        function obj = simDriver(timestep, truth_mat_file_path, main_filter, ...
                                  backup_filter, parachute_deployment_logic, output_dir_path)

            %SIM_DRIVER Construct an instance of this class
            %   Detailed explanation goes here

            % Timestep in seconds
            obj.timestep = timestep;

            % Filters
            obj.main_filter = main_filter;
            obj.backup_filter = backup_filter;
            obj.parachute_deployment_logic = parachute_deployment_logic;
            obj.output_dir_path = output_dir_path;

            % Truth Documents (sourced from Flight D vehicle trajectory sims)
            obj.truth_file_path = truth_mat_file_path;
            obj.truth = load(truth_mat_file_path);
            obj.truth_hash = hash_mat_file(obj.truth);

            % Sensors Container 
            obj.sensors = containers.Map( ...
                'KeyType', 'char', ...
                'ValueType', 'any');        
        end

        function addSensor(obj, sensorType, sensor)
            %{
                Example Usage:
    
                imu = ASM330LHBG1(1000, accelParams, gyroParams);
                obj.addSensor('IMU', imu);
            %}

            key = char(string(sensorType));

            stride = (1 / sensor.sampleRate) / obj.timestep;
            strideRounded = round(stride);

            if strideRounded < 1 || abs(stride - strideRounded) > 1e-8
                error("sim_driver:InvalidSampleRate", ...
                    "Sensor sample period must be an integer multiple " + ...
                    "of the simulation timestep.");
            end

            obj.sensors(key) = sensor;
            obj.sensor_gather_intervals(key) = strideRounded;
            obj.sensor_measurement_types = [obj.sensor_measurement_types; sensor.sensor_category];       
        end

        function runSimulation(obj)
            
            simData = obj.truth.simData;        % Get the SimData struct
            simulationTime = simData.time(:);   % Get the time array
            numSteps = numel(simulationTime);   % Get the number of timesteps
            sensorKeys = keys(obj.sensors);     % 

            %% Timestep and Truth Delta T Must Be the Same
            truth_dt = simData.time(2) - simData.time(1);   % Determine the timestep of the truth/SimData struct

            if (truth_dt ~= single(obj.timestep))
                error("Simulator timestep is not the same as FlightD timestep")
            end

            %% Sensor Gather Intervals
            totalMeasurements = 0;      % Total Number of Measurements from all sensors 

            for sensorIndex = 1:numel(sensorKeys)
                sensor = obj.sensors(sensorKeys{sensorIndex}); % Get a sensor
                totalMeasurements = totalMeasurements + sensor.num_measurements; % Add its measurements to the total measurements
            end

            currentMeasurement = zeros(totalMeasurements, 1); % 
            
            %% Simulation
            for stepIndex = 1:numSteps
                currentTime = simulationTime(stepIndex); % Get the current time within the simulation
                currMeasIndex = 1;                       % Initialize the current measurement index to the start of currentMeasurement

                % Loop through the sensors
                for sensorIndex = 1:numel(sensorKeys)
                    key = sensorKeys{sensorIndex};              % Get the name/key of the sensor
                    sensor = obj.sensors(key);                  % Get the actual sensor object
                    stride = obj.sensor_gather_intervals(key);  % The stride is the number of timesteps required before a sensor measurement is valid
                    sensorNumMeas = sensor.num_measurements;    % Number of measurements this sensor outputs
                    
                    if mod(stepIndex - 1, stride) == 0
                        % We have a measurement from this sensor. Go to the following if statements
                        % to determine what kind of sensor this is.

                        if (sensor.sensor_category == measurement_type.IMU)
                            %% IMU Measurment (NOTE: SHOULD OCCUR AT EVERY TIME STEP)

                            trueAccel = simData.acc(:,stepIndex); % Get the true acceleration from the simData
                            trueGyro = simData.wdot(:,stepIndex); % Get the true angular velocity (wdot) from the simData
                            currentMeasurement(currMeasIndex : currMeasIndex + sensorNumMeas) = sensor.measure(trueAccel, trueGyro); % Get the measured accel and wdot and put it in the currentMeasurement

                        elseif (sensor.sensor_category == measurement_type.BARO)
                            %% Barometer Measurement

                            wgs84_altitude = simData.lla0(3,1) + simData.pos(3,simIndex);  % Get the altitude by adding the inital altitude with respect to the WGS84 ellipsoid to the current altitude
                            [~, ~, truePressure, ~] = atmosisa(wgs84_altitude);            % Use the International Standard Atmosphere to turn the altitude into a pressure
                            currentMeasurement(currMeasIndex : currMeasIndex + sensorNumMeas) = sensor.measure(truePressure);   % Get the actual pressure measured by the barometer and put it into currentMeasurment

                        elseif (sensor.sensor_category == measurement_type.MAG)
                            %% Magnetometer Measurement


                        elseif (sensor.sensor_category == measurement_type.GPS)
                            %% GPS Measurement

                        end
                    else
                            %% No Measurement 

                        % We don't have a measurement from this sensor and therefore fill its array
                        % with NaN (aka Not a Number)
                        currentMeasurement(currMeasIndex : currMeasIndex + sensorNumMeas) = NaN(sensorNumMeas, 1);
                    end
                end
                
                %% Running the Filter

                mainFilter.runTimestep(currentMeasurement, obj.sensor_measurement_types, dt)

            end

            %% Simulation Has Ended

        end
    end
end
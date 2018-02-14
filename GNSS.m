classdef GNSS < handle
   properties(Constant)
       % Constants
        deg_to_rad = 0.01745329252; % Degrees to radians conversion factor
        rad_to_deg = 1/0.01745329252; % Radians to degrees conversion factor
        c = 299792458; % Speed of light in m/s
        omega_ie = 7.292115E-5;  % Earth rotation rate in rad/s
        Omega_ie = Skew_symmetric([0,0,7.292115E-5]);
        R_0 = 6378137; %WGS84 Equatorial radius in meters
        e = 0.0818191908425; %WGS84 eccentricity
   end
   properties
       pseudo_range_data; % Raw data for pseudo-range
       pseudo_range_rate_data;     
       num_satellite; % number of satellite
       num_time; % number of time
       time_data; % Measurement time
       satellite_number; % satellite index
       time_index = 1; % current time index
       %satellite_index = 1; % current satellite
       pseudo_ranges_est; % Esitimated pseudo-ranges
       pseudo_range_rates_est;
       satellite_positions; % Satellite Position, 3 by 8
       satellite_velocities; % Satellite velocities, 3 by 8 
       LOS; % Line-of-Sight unit vector
       
       % State vector
       user_pos_ned;
       user_pos_ecef;
       user_vel_ned = 0;
       user_vel_ecef = 0;
       clock_offset = 0;
       clock_drift= 0;
   end
   
   methods
       function obj = GNSS()
            pos_data = csvread('Pseudo_ranges.csv');
            speed_data = csvread('Pseudo_range_rates.csv');
            obj.time_data = pos_data(2:size(pos_data,1),1);
            obj.satellite_number = pos_data(1, 2:end);
            obj.pseudo_range_data = pos_data(2:end, 2:end);
            obj.pseudo_range_rate_data = speed_data(2:end, 2:end);
            obj.num_satellite = length(obj.satellite_number);
            obj.num_time = length(obj.time_data);
       end
       
       % Compute satellite position and velocity
       function calculate_Satellite_Position_Vel(obj)
            obj.satellite_positions = zeros(3, obj.num_satellite);
            obj.satellite_velocities = zeros(3, obj.num_satellite);
            t = obj.time_data(obj.time_index);
            for n = 1 : obj.num_satellite
                [obj.satellite_positions(:, n), obj.satellite_velocities(:, n)] = Satellite_position_and_velocity(t, obj.satellite_number(n));
            end
       end
       
       % Estimate the pseudo-range given the estimated user position
       function est_Pseudo_range(obj)
           obj.pseudo_ranges_est = zeros(obj.num_satellite, 1);
           for n = 1 : obj.num_satellite
                curr_p = 0;
                prev_p = ones(3,1);
                while norm(curr_p - prev_p) >= 0.0005
                    prev_p = curr_p;
                    C = eye(3);
                    C(1,2) = obj.omega_ie * curr_p / obj.c;
                    C(2,1) = -obj.omega_ie * curr_p / obj.c;
                    curr_p = sqrt((C * obj.satellite_positions(:,n) - obj.user_pos_ecef)' * (C * obj.satellite_positions(:,n) - obj.user_pos_ecef));
                end
                obj.pseudo_ranges_est(n) = curr_p;
           end
       end
       
       
       function est_Pseudo_range_rate(obj)
           obj.pseudo_range_rates_est=zeros(obj.num_satellite, 1);
           for n = 1 : obj.num_satellite
                C = eye(3);
                C(1,2) = obj.omega_ie * obj.pseudo_ranges_est(n) / obj.c;
                C(2,1) = -obj.omega_ie * obj.pseudo_ranges_est(n) / obj.c;
                u_aj= obj.LOS(:,n)';
                r_ej = obj.satellite_positions(:, n);
                r_ea = obj.user_pos_ecef;
                v_ej = obj.satellite_velocities(:, n);
                v_ea = obj.user_vel_ecef;
                Omega = obj.Omega_ie;
                obj.pseudo_range_rates_est(n)=u_aj*(C*(v_ej+Omega*r_ej)-(v_ea+Omega*r_ea));
           end
       end
       
       
       
       
       % Compute Line-of-Sight Unit Vector
       function compute_LineofSight(obj)
            obj.LOS = zeros(3, obj.num_satellite); 
            for n = 1 : obj.num_satellite
                 obj.LOS(:, n) = (obj.satellite_positions(:,n) - obj.user_pos_ecef) / obj.pseudo_ranges_est(n);
            end
       end
       
       % Least square estimation
       function least_square_pos(obj)
            line_of_sight = obj.LOS;
            x_pos = [obj.user_pos_ecef; obj.clock_offset];
            %x_vel = [obj.user_vel_ecef; obj.clock_drift];
            Hg = [-line_of_sight', ones(obj.num_satellite, 1)];
            pseudo_ranges = obj.pseudo_range_data(obj.time_index, :);
            %pseudo_range_rates = obj.pseudo_range_rate_data(obj.time_index, :);
            deltaZ_pos = pseudo_ranges' - obj.pseudo_ranges_est - obj.clock_offset;
            %deltaZ_vel = pseudo_range_rates' - obj.pseudo_range_rates_est
            %- obj.clock_drift
            x_plus_pos = x_pos + (Hg' * Hg) \ Hg' * deltaZ_pos;
            obj.user_pos_ecef = x_plus_pos(1:3);
            %obj.user_vel_ecef = x_plus_vel(1:3);
            obj.clock_offset= x_plus_pos(4);
       end
       
       function least_square_vel(obj)
           line_of_sight = obj.LOS;
            x_vel = [obj.user_vel_ecef; obj.clock_drift];
            Hg = [-line_of_sight', ones(obj.num_satellite, 1)];
            %pseudo_ranges = obj.pseudo_range_data(obj.time_index, :);
            pseudo_range_rates = obj.pseudo_range_rate_data(obj.time_index, :);
            %deltaZ_pos = pseudo_ranges' - obj.pseudo_ranges_est - obj.clock_offset;
            deltaZ_vel = pseudo_range_rates' - obj.pseudo_range_rates_est- obj.clock_drift;
            x_plus_vel = x_vel + (Hg' * Hg) \ Hg' * deltaZ_vel;
            %obj.user_pos_ecef = x_plus_pos(1:3);
            obj.user_vel_ecef = x_plus_vel(1:3);
            obj.clock_drift= x_plus_vel(4);
       end
       
       
       
   end      
end
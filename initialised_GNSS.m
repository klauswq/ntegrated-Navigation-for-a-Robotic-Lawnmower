function [x_est,P_matrix] = initialised_GNSS
    Define_Constants;

    pos_data = csvread('Pseudo_ranges.csv');
    speed_data = csvread('Pseudo_range_rates.csv');
    time_data = pos_data(2:size(pos_data,1),1);
    satellite_number = pos_data(1, 2:end);
    pseudo_range_data = pos_data(2:end, 2:end);
    pseudo_range_rate_data = speed_data(2:end, 2:end);
    num_satellite = length(satellite_number);

    clock_offset = 0;
    clock_drift= 0;

    user_pos_ecef = zeros(3,1);
    user_vel_ecef = zeros(3,1);
    prev_pos = ones(3,1);

    while norm(prev_pos(1:3) - user_pos_ecef) > 0.001 
        prev_pos = user_pos_ecef;    

        satellite_positions = zeros(3, num_satellite);
        satellite_velocities = zeros(3, num_satellite);
        t = time_data(1);
        for n = 1 : num_satellite
            [satellite_positions(:, n), satellite_velocities(:, n)] = Satellite_position_and_velocity(t, satellite_number(n));
        end

        pseudo_ranges_est = zeros(num_satellite, 1);
        for n = 1 : num_satellite
            curr_p = 0;
            prev_p = ones(3,1);
            while norm(curr_p - prev_p) >= 0.0005
                prev_p = curr_p;
                C = eye(3);
                C(1,2) = omega_ie * curr_p / c;
                C(2,1) = -omega_ie * curr_p / c;
                curr_p = sqrt((C * satellite_positions(:,n) - user_pos_ecef)' * (C * satellite_positions(:,n) - user_pos_ecef));
            end
            pseudo_ranges_est(n) = curr_p;
        end

        % Compute Line-of-Sight Unit Vector
        LOS = zeros(3, num_satellite); 
        for n = 1 : num_satellite
            LOS(:, n) = (satellite_positions(:,n) - user_pos_ecef) / pseudo_ranges_est(n);
        end

        % Least square estimation
        line_of_sight = LOS;
        x_pos = [user_pos_ecef; clock_offset];

        Hg = [-line_of_sight', ones(num_satellite, 1)];
        pseudo_ranges = pseudo_range_data(1, :);
        deltaZ_pos = pseudo_ranges' - pseudo_ranges_est - clock_offset;

        x_plus_pos = x_pos + (Hg' * Hg) \ Hg' * deltaZ_pos;
        user_pos_ecef = x_plus_pos(1:3);

        clock_offset= x_plus_pos(4);
    end

    pseudo_range_rates_est=zeros(num_satellite, 1);
    for n = 1 : num_satellite
        C = eye(3);
        C(1,2) = omega_ie * pseudo_ranges_est(n) / c;
        C(2,1) = -omega_ie * pseudo_ranges_est(n) / c;
        u_aj= LOS(:,n)';
        r_ej = satellite_positions(:, n);
        r_ea = user_pos_ecef;
        v_ej = satellite_velocities(:, n);
        v_ea = user_vel_ecef;
        pseudo_range_rates_est(n)=u_aj*(C*(v_ej+Omega_ie*r_ej)-(v_ea+Omega_ie*r_ea));
    end

    line_of_sight = LOS;
    x_vel = [user_vel_ecef; clock_drift];
    Hg = [-line_of_sight', ones(num_satellite, 1)];

    pseudo_range_rates = pseudo_range_rate_data(1, :);
    deltaZ_vel = pseudo_range_rates' - pseudo_range_rates_est- clock_drift;
    x_plus_vel = x_vel + (Hg' * Hg) \ Hg' * deltaZ_vel;

    user_vel_ecef = x_plus_vel(1:3);
    clock_drift= x_plus_vel(4);

    % the position of the user is fixed, now we will compute the velocity and
    % clock offset and clock drift. initially, we assume that the first
    % velocities of the user to be zeros
    x_est = [user_pos_ecef(1);user_pos_ecef(2);user_pos_ecef(3);...
             user_vel_ecef(1); user_vel_ecef(2); user_vel_ecef(3);...
             clock_offset; clock_drift];

    % Initialise error covariance matrix
    P_matrix =  zeros(8);
    P_matrix(1,1) = 100;
    P_matrix(2,2) = 100;
    P_matrix(3,3) = 100;
    P_matrix(4,4) = 0.01;
    P_matrix(5,5) = 0.01;
    P_matrix(6,6) = 0.01;
    P_matrix(7,7) = clock_offset;
    P_matrix(8,8) = clock_drift;
end



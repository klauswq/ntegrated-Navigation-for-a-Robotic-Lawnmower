function result = Integration

    Define_Constants;
    pseudo_range = load('Pseudo_ranges.csv');
    Time = pseudo_range(2:end,1);

    % *** Step 1: Get GNNS positions
    GNSSResult = ComputeGNSS;
    [num_epoch,n] = size(GNSSResult);
    GNSSResult(1,3:4) =zeros(1,2) ;
    gn_position = GNSSResult(:,1:2);
    gn_velocity = GNSSResult(:,3:4);
    gn_heading =  GNSSResult(:,5); 

    % *** Step 2: Get dead reckoning solutions
    DKResult = ComputeDeadreckoning;
    dk_position = DKResult(:,1:2);
    dk_velocity = DKResult(:,3:4);
    dk_heading =  DKResult(:,5);

    pos = (dk_position'*gn_position)./num_epoch;
    vel = (dk_velocity'*gn_velocity)./num_epoch;
    head = (dk_heading'*gn_heading)./num_epoch * deg_to_rad;

    meas_noise_matrix = [pos, zeros(2), zeros(2,1); 
                         zeros(2),vel,  zeros(2,1);
                         zeros(1,2), zeros(1,2),head];

    % *** Step 3: Kalman filter  
    [x_est,P_matrix] = Initialise_Integration_KF;

    % results:
    % the first column is Time in seconds
    % the second column is Geodetic latitude in degrees
    % the third column is Geodetic longitude in degrees
    % the fourth column is North velocity in metres per second
    % the fifth column is East velocity in metres per second
    % the sixth column is Heading in degrees

    results = zeros(num_epoch, 6);
    calibration_values = zeros(num_epoch, 5);
    for i=1:num_epoch
        inte_dk = DKResult(i,:);
        inte_GNSS = GNSSResult(i,:);

        % get the latitude in this epoch, use it as the input parameter to
        % integration kalman filter
        latitude =  inte_GNSS(1)*deg_to_rad;

        % the previous latitude data will be used when calculate the transition
        % matrix
        if i>1
            pre_latitude = GNSSResult(i-1,1)*deg_to_rad;
        else
            pre_latitude = latitude;
        end

        [newState, new_error_cov] = kalmanFilter(x_est, P_matrix,inte_dk,inte_GNSS,...
                                                 meas_noise_matrix,latitude, pre_latitude);
        calibration_values(i,:) = newState(1:5)';
        calibration_values(i,5) = calibration_values(i,5)*rad_to_deg;

        %the first column of the results is the time in seconds 
        results(i,1) = Time(i);

        %the output of the integration kalman filter will be used to calibrate the DR result 
        results(i,2:6) = inte_dk - calibration_values(i,:);

        x_est = newState;
        P_matrix = new_error_cov;
    end

    figure;
    result = results ;
    geoshow(result(:,2),result(:,3));
    hold off;
    figure;
    plot(Time, result(:,4),'r-o');
    hold on;
    plot(Time, result(:,5),'b-x');
    dlmwrite('map.csv',[result(:,2),result(:,3)],'precision','%.6f');
    dlmwrite('Output_Profile.csv',[result(:,1),result(:,2),result(:,3),result(:,4),result(:,5),result(:,6)],'precision','%.6f');

end

function [newState, new_error_cov] = kalmanFilter(x_est, P_matrix,odometry_sol,GNSSValue, meas_noise_matrix,latitude, pre_latitude)
    Define_Constants ;
    prog_interval = 0.5 ;
    [R_N,R_E] = Radii_of_curvature(latitude);
    %  *** Step 1: Compute transition matrix

    v = [prog_interval/R_N , prog_interval/(R_E*cos(pre_latitude))];
    v_error = diag(v);

    % x = [pos_lat, pos_lon, vel_N, vel_E, heading, gyro_angular_vel, gyro_bias]'
    transMatrix = [eye(2), v_error, zeros(2,1)  , zeros(2,1), zeros(2,1) ; 
                   zeros(2), eye(2),zeros(2,1), zeros(2,1), zeros(2,1)   ;    
                   zeros(1,2), zeros(1,2), 1, prog_interval, 0 ;
                   zeros(1,2), zeros(1,2), 0, 1, prog_interval;
                   zeros(1,2), zeros(1,2), 0, 0, 1] ;

    % Step 2: Compute system noise covariance matrix 

    % space_err_std = 1;
    % Residual_ionosphere_err_std = 2;
    % Residual_troposphere_err_std = 0.2;
    % Code_tracking_multipath_err_std = 2;
    % Range_rate_tracking_multipath_err_std = 0.02;
    % Receiver_clock_drift_std = 200;
    % Wheel_Scale_factor_err_std = 0.03;
    % Wheel_noise_std = 0.05;
    Gyro_bias_std = 1 * deg_to_rad;
    Gyro_scale_factor_err_std = 0.01;
    Gyro_cross_coupling_err_std = 0.001;
    Gyro_noise_std = 0.0001;

    Gyro_bias_std = 1 * deg_to_rad;

    S_heading = (2*deg_to_rad)^2;
    %S_gyro = (3*10^4)^2/prog_interval ;
    S_bgd = ((Gyro_bias_std/prog_interval)^2)/prog_interval;

    Wheel_Scale_factor_err_std = 0.03;
    Wheel_noise_std = 0.05;              
    wheel_noise = Wheel_Scale_factor_err_std + Wheel_noise_std;
    gyro_noise = Gyro_bias_std + Gyro_scale_factor_err_std + ...
                 Gyro_cross_coupling_err_std + Gyro_noise_std;

    Q_sys_noise = [zeros(2), zeros(2), zeros(2,1)  , zeros(2,1), zeros(2,1);...
                   zeros(2), wheel_noise*prog_interval*eye(2), zeros(2,1), zeros(2,1), zeros(2,1); ...
                   zeros(1,2), zeros(1,2), S_heading, 0,0; ...
                   zeros(1,2), zeros(1,2), 0, gyro_noise*prog_interval, 0;
                   zeros(1,2), zeros(1,2),0 , 0, S_bgd*prog_interval];

    % *** Step 3: Propagate the state 
    new_x_est = transMatrix * x_est;

    % *** Step 4: Propagate the covariance
    error_P_matrix =  transMatrix * P_matrix * transMatrix.' + Q_sys_noise;

    % *** Step 5: Compute the measurement matrix 
    meas_matrix = [-eye(2), zeros(2),  zeros(2,1), zeros(2,1), zeros(2,1); ...
                    zeros(2), -eye(2), zeros(2,1)  , zeros(2,1), zeros(2,1); 
                    zeros(1,2), zeros(1,2), -1, 0, 0];

    % *** Step 6: Compute the measurement noise covariance matrix 
    %we computed it before with the variance between measurement 
    noise_cov = meas_noise_matrix;

    % *** Step 7: Compute Kalman gain matrix 
    K_gain = error_P_matrix*meas_matrix.'*(meas_matrix*error_P_matrix*meas_matrix.'+noise_cov)^-1;

    % *** Step 8: Compute the measurement innovation
    meas_Inv = [GNSSValue(1:2)'-odometry_sol(1:2)' + new_x_est(1:2);
                       GNSSValue(3:4)'-odometry_sol(3:4)' + new_x_est(3:4);
                       GNSSValue(5)*deg_to_rad-odometry_sol(5)*deg_to_rad + new_x_est(5)];

    newState = new_x_est + K_gain * meas_Inv ;
    new_error_cov = (eye(7)- K_gain * meas_matrix) * error_P_matrix;
end



function GNSS_Result=ComputeGNSS
    Define_Constants;
    pseudo_range=csvread('Pseudo_ranges.csv');
    pseudo_range_rate=csvread('Pseudo_range_rates.csv');
    satellites = pseudo_range(1, 2:end);
    satellites_inTotal=size(satellites,2);
    satellite_Pos_Vel = zeros(3, satellites_inTotal,2);
    num_epoch = size(pseudo_range,1)-1;
    time = pseudo_range(2:end,1);

    %store the poisition and velocity data of all the satellites at epoch 0
    for i=1:satellites_inTotal
        this_satellite=satellites(i);
        [position, velocity]=Satellite_position_and_velocity(0,this_satellite);
        satellite_Pos_Vel(:,i,1)=position;
        satellite_Pos_Vel(:,i,2)=velocity;
    end
    
    %initialisation 
    [x_est, P_matrix]=initialised_GNSS;
  
    GNSS_Result=zeros(num_epoch,5);
    
    % store epoch 1 data into GNSS_Result
    [L_b, lambda, ~, v_eb_n]=pv_ECEF_to_NED(x_est(1:3),x_est(4:6));
    
    lat_predic = L_b * rad_to_deg;   
    long_predic = lambda * rad_to_deg;   
    GNSS_Result(1,1:4) = [lat_predic, long_predic, v_eb_n(1:2)'];
    
    for i = 2:num_epoch
        measured_ranges = pseudo_range(i+1, 2:end);      
        measured_range_rates = pseudo_range_rate(i+1, 2:end);      
        time_index=pseudo_range(i+1,1);
        
        % estimate the new state and cov matrix using a kalman filter        
        [new_est_state,new_error_cov,outlier] = kalman_filter(measured_ranges, measured_range_rates, x_est, P_matrix, time_index, satellites);
    
        % take out the outliers and recompute the position        
        if ~isempty(outlier)
            [new_range, new_range_rate, new_satellites] = outlier_removed(outlier,measured_ranges, measured_range_rates, satellites);
            [new_est_state,new_error_cov,~]=kalman_filter(new_range, new_range_rate, x_est, P_matrix, time_index, new_satellites);
        end
        
        % transfer ECEF_user pos and vel into NED coordinator
        [L_b,lambda_b,~,v_eb_n] = pv_ECEF_to_NED(new_est_state(1:3),new_est_state(4:6)) ;
        
        latPredic = L_b*rad_to_deg ;
    	longPredic = lambda_b*rad_to_deg ;
        GNSS_Result(i,1:4) = [latPredic, longPredic, v_eb_n(1:2)'];
        
        % compute the heading
        heading = atan(GNSS_Result(i,4)/ GNSS_Result(i,3)).*rad_to_deg ;
        if GNSS_Result(i,3) < 0
            heading = heading + 180;
        end
        GNSS_Result(i,5) = heading;
        
        x_est = new_est_state;
        P_matrix = new_error_cov;       
    end
    show = 1;
    if (show)
    % *** Show the position
    figure
        geoshow(GNSS_Result(:,1),GNSS_Result(:,2));
        title('GNSS with outliers');
        xlabel('longitude');
        ylabel('lattidue');

    % *** Show the velocity
    figure
    plot(time, GNSS_Result(:,3),'r-o');
    hold on;
    plot(time, GNSS_Result(:,4),'b-x');   
    end

    dlmwrite('GNSS_Output.csv',GNSS_Result,'precision','%.6f');    
end
       
    function [new_est_state,new_error_cov,outlier] = kalman_filter(ranges_m, range_rates_m, x_est, P_matrix, time_index, satellites)
        Define_Constants;
        
        %step 1 Compute the transition matrix 
        tau=0.5;

        transition_matrix=[eye(3,3), tau*eye(3,3), zeros(3,1), zeros(3,1);...
                           zeros(3,3), eye(3,3), zeros(3,1), zeros(3,1);... 
                           zeros(1,3), zeros(1,3), 1, tau;...
                           zeros(1,3), zeros(1,3), 0, 1];
               
        %step 2 Compute the system noise covariance matrix
        
        S_ae = 1; % acceleration PSD
        S_phi = 0.01;  %Clock phase PSD: 0.01 
        S_f=0.04; %Clock frequency PSD: 0.04
        
        System_Noise_Cov=[1/3*S_ae*tau^3*eye(3), 1/2*S_ae*tau^2*eye(3), zeros(3,1), zeros(3,1) ;...
                          1/2*S_ae*tau^2*eye(3), S_ae*tau*eye(3), zeros(3,1), zeros(3, 1); ...
                          zeros(1,3), zeros(1,3), S_phi*tau+1/3*tau^3*S_f, 1/2*tau^2*S_f ; ...
                          zeros(1,3), zeros(1,3),  1/2*tau^2*S_f, S_f*tau] 
                      
        
        
        %step 3 propagate the state estimation
        
        x_est_updated = transition_matrix * x_est;
        
        %step 4 propogate the error covariance estimation
        
        error_cov = transition_matrix * P_matrix * transition_matrix' + System_Noise_Cov;
        
        %step 5 compute the measurement matrix
        
        num_satellites=size(satellites,2);
        
        satellite_positions = zeros(3, num_satellites);
        
        satellite_velocities = zeros(3, num_satellites);
        
        r_aj=zeros(1,num_satellites);
        
        for i = 1: num_satellites
            this_satellite=satellites(1,i);
            [satellite_positions(:,i), satellite_velocities(:,i)] = Satellite_position_and_velocity(time_index, this_satellite);
            p = 0;
            prev_p = ones(3,1);
            while norm(p - prev_p) >= 0.001
                prev_p = p;
                C = eye(3);
                C(1,2) = omega_ie * p / c;
                C(2,1) = -omega_ie * p / c;
                quad = (C * satellite_positions(:,i) - x_est_updated(1:3,:))' * (C * satellite_positions(:,i) - x_est_updated(1:3,:));
                p = sqrt(quad);
            end
            r_aj(1,i) = p;
        end
        
        %compute the unite of line_of_sight
        
        LineOfSight = (satellite_positions - repmat(x_est_updated(1:3,:), [1,num_satellites]))./repmat(r_aj,[3,1]);
        
        %measurement matrix
        
        measured_matrix=[-LineOfSight.' zeros(num_satellites, 3) ones(num_satellites,1) zeros(num_satellites,1) ; ...
                         zeros(num_satellites, 3) -LineOfSight.' zeros(num_satellites,1) ones(num_satellites,1)] ;
                     
        %step 6 Compute the measurement noise covariance matrix
        pseudo_range_error = 10;
        pseudo_range_rate_error = 0.05;
        R_k = [pseudo_range_error^2 * eye(num_satellites,num_satellites),zeros(num_satellites,num_satellites);
               zeros(num_satellites,num_satellites),pseudo_range_rate_error^2 * eye(num_satellites,num_satellites)];

        
        %step 7 Compute the Kalman gain matrix
        
        Kalman_gain = error_cov * measured_matrix.'*(measured_matrix * error_cov * measured_matrix.' + R_k)^(-1)
        
        r_aj_dot = zeros(1, num_satellites);
        for i=1:num_satellites
            % the measured pos and vel for satellite i
            satPos = satellite_positions(:,i);
            satVel = satellite_velocities(:,i);
            % propogated state 
            user_pos_predic = x_est_updated(1:3,:);
            user_vel_predic = x_est_updated(4:6,:);
            % predicted pseudo range
            range = r_aj(1,i);
            
            % predict the range rate
            C = [1, omega_ie*range/c, 0 ; -omega_ie*range/c 1 0 ; 0 0 1] ;
            r_aj_dot(1,i) = LineOfSight(:,i).' * (C*(satVel+Omega_ie*satPos)-(user_vel_predic+Omega_ie*user_pos_predic));
           
        end
        
        % step 8 Formulate the measurement innovation vector
        
        p_a_j = ranges_m';
        p_a_j_dot = range_rates_m';
        p_c_a = x_est_updated(7,:);
        p_c_a_dot = x_est_updated(8,:);
        
        measurement_innov = [p_a_j - r_aj' - p_c_a ;p_a_j_dot - r_aj_dot' - p_c_a_dot];
        
        % step 9 Compute the new state
        
        new_est_state = x_est_updated + Kalman_gain * measurement_innov;
        
        % step 10 Compute the new covariance matrix
        
        new_error_cov = (eye(8) - Kalman_gain * measured_matrix) * error_cov;
        
        % Outliers detection
        outlier = outlierDetector(measured_matrix,measurement_innov,LineOfSight,num_satellites);       
    end
    
    function outlier = outlierDetector(m_matrix,m_innov,LOS,num_sat)
        % Compute the residuals vector
        
        m_matrix_pos=[m_matrix(1:num_sat, 1:3), ones(num_sat,1)];
        
        m_innov_pos = m_innov(1:num_sat,:);
        
        residual = (m_matrix_pos*(m_matrix_pos'*m_matrix_pos)^(-1)*m_matrix_pos'-eye(num_sat))*m_innov_pos;
        
        % Compute the residuals covariance matrix
        
        measErrorSD = repmat(2.2, [1,num_sat ]) ; %2.2=2.0+0.2
        measErrorVariance = (measErrorSD./LOS(3,:)).^2 ;
        
        C_v=(eye(num_sat)-m_matrix_pos*(m_matrix_pos'*m_matrix_pos)^(-1)*m_matrix_pos');
        
        % Compute the diagonal element of C_v
        
        diag_C_v=diag(C_v)'.*measErrorVariance;
        
        % Compute the normalized residual
        
        norm_residual = residual'./sqrt(diag_C_v);
        
        threshold=6;
        
        T = abs(norm_residual)>(diag_C_v).*threshold;
        
        if sum(T) > 0
            outlier = T ;
        else
            outlier = [] ;
        end
    end


    function [new_range, new_range_rate, new_satellites] = outlier_removed(outlier,measured_ranges, measured_range_rates, satellites)

        number_Sat = size(satellites,2) ;
        new_range = measured_ranges ;
        new_range_rate = measured_range_rates ;
        new_satellites = satellites ;
        
        % if a satellite is an outlier, we ignore its value
        index = 0 ;
        for i=1:number_Sat
            if outlier(i)==1
                index = index+1 ;
                new_range(i-index) = [] ;
                new_range_rate(i-index) = [] ;
                new_satellites(i-index) = [] ;
            end
        end       
    end
        

    
    
    
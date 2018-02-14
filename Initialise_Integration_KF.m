function [x_est,P_matrix] = Initialise_Integration_KF
Define_Constants ;

% Initialise state estimates
x_est = zeros(7,1);

% Initialise error covariance matrix
P_matrix =  zeros(7);
P_matrix(1:2,1:2) = eye(2) * 0.01^2; % position error to 0.01 rad
P_matrix(3:4,3:4) = eye(2) * 0.01^2; % velocity error to 0.1m/s 
P_matrix(5,5) = (2*deg_to_rad)^2; % initialize the heading to 2 degree
P_matrix(6,6) = (1*deg_to_rad)^2 ; % gyro angular velocity 1deg/s
P_matrix(7,7) = (1*deg_to_rad)^2 ; % gyro bias 1 degree per second 
end
function result = ComputeDeadreckoning
    clear
    Define_Constants
    Dead_reckoning = load('Dead_reckoning.csv');
    [num_epoch,n] = size(Dead_reckoning);
    time = Dead_reckoning(:,1);
    % *** columns 2 to 5 contain the wheel-speed measurements
    wheel_speed = Dead_reckoning(:,2:5);
    % the rear wheels are the driving wheels
    Rear_l = wheel_speed(:,3);
    Rear_r = wheel_speed(:,4);
    a_speed = (Rear_l + Rear_r)/2;

    %heading = Dead_reckoning(:,7)*deg_to_rad
    heading = newheading;

    % *** Step 1: the average velocity between each epochs 
    new_vel = zeros(num_epoch,2);
    for i = 2:num_epoch
        new_vel(i-1,:) = 1/2 * a_speed(i,1) * [cos(heading(i,1))+cos(heading(i-1,1)) sin(heading(i,1))+sin(heading(i-1,1))];
    end

    % *** Step 2: compute position from counterparts at last epoch
    position = zeros(num_epoch,2);
    latitude = 51.5092544; % from GNSS reault
    longitude = -0.1610454;
    height = -0.0;
    position(1,1) = latitude*deg_to_rad;
    position(1,2) = longitude*deg_to_rad;

    % R_N is the meridian radius of curvature
    % R_E is the transversr radius of curvature
    [R_N,R_E] = Radii_of_curvature(latitude*deg_to_rad);
    for i = 2:num_epoch
            position(i,1) = position(i-1,1) + (new_vel(i-1,1) * 0.5)/(R_N+height); %north
            position(i,2) = position(i-1,2) + (new_vel(i-1,2) * 0.5)/((R_E+height)*cos(position(i,1))); %east
    end

    % *** Step 3: compute the instantaneouse DR velocity at each epoch
    position = position*rad_to_deg;
    V_N = zeros(num_epoch,1);
    V_E = zeros(num_epoch,1);
    V_N(1,1) = a_speed(1,1)*cos(heading(1,1));
    V_E(1,1) = a_speed(1,1)*sin(heading(1,1));
    for i = 2:num_epoch
        V_N(i,1) = 2*new_vel(i-1,1)-V_N(i-1,1);
        V_E(i,1) = 2*new_vel(i-1,2)-V_E(i-1,1);
    end

    show = 1;
    if (show)
    % *** Show the position
    figure
        geoshow(position(:,1),position(:,2));
        title('Deadreckoning');
        xlabel('longitude');
        ylabel('lattidue');

    % *** Show the velocity
    figure
    plot(time, new_vel(:,1),'r-o');
    hold on;
    plot(time, new_vel(:,2),'b-x');


    % *** Show the heading
    figure
        plot(time,heading);
        title('Deadreckoning heading');
        xlabel('time(s)');
        ylabel('heading');
    end
    result = [position,new_vel,heading]
    % *** Gen the output file
    dlmwrite('Deadreckoning_Output.csv',[position,new_vel(:,1),new_vel(:,2),heading*rad_to_deg],'precision','%.6f');
end
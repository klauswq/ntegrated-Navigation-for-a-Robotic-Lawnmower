%% Combine Gyroscope-derived heading and Magnetic heading 
function new_heading = newheading   
    Define_Constants;
    Dead_reckoning = load('Dead_reckoning.csv');

    gyroscope = Dead_reckoning(:,6);
    heading = Dead_reckoning(:,7) * deg_to_rad;
    time = Dead_reckoning(:,1);
    [num_epoch,n] = size(Dead_reckoning);    

    new_heading = zeros(num_epoch,1);
    new_heading(1,:) = heading(1,:);
   
    K = 0.5;   
    for i = 2:num_epoch
        new_heading(i,:) = K * heading(i,:) + (1-K) * (new_heading(i-1,:)+ gyroscope(i,:) * (time(i)-time(i-1)));
    end
 end 

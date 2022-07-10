%% Atitude Estimation Complementary Filter 
% Acc/Mag CF gyro (Acc/Mag)
% inputs:   Accelerometer Measurements  Acc:  n x 3 matric
%           Magnetometer Measurements   Mag:  n x 3 matric
%           Reference Quaternion        Quat: n x 4 matric
%           Sample Rate Frequency       fs:   1 x 1 scaler (Hz)
function [RMSE_CF] = CFgyroCF(Acc,Gyro,Mag,Quat,fs)
 for i = 1:length(Acc)
    acc_x(i,1) = Acc(i,1)/norm(Acc(i,:));
    acc_y(i,1) = Acc(i,2)/norm(Acc(i,:));
    acc_z(i,1) = Acc(i,3)/norm(Acc(i,:));
    mag_x(i,1) = Mag(i,1)/norm(Mag(i,:));
    mag_y(i,1) = Mag(i,2)/norm(Mag(i,:));
    mag_z(i,1) = Mag(i,3)/norm(Mag(i,:));
    if (norm(Gyro(i,:))) == 0
        gyro_x(i,1) =   Gyro(i,1);
        gyro_y(i,1) =   Gyro(i,2);
        gyro_z(i,1) =   Gyro(i,3);
    else
        gyro_x(i,1) =   Gyro(i,1)/norm(Gyro(i,:));
        gyro_y(i,1) =   Gyro(i,2)/norm(Gyro(i,:));
        gyro_z(i,1) =   Gyro(i,3)/norm(Gyro(i,:));
    end
end   
% calculate the time
    dt      =   1/fs;
    t       =   0:dt:(length(Acc) - 1)*dt;
%% Calculate Attitude Reference
ref_eul = quat2eul(Quat);
for i = 1:length(ref_eul)
    phi_ref(i,1) = ref_eul(i,3);
    theta_ref(i,1)=ref_eul(i,2);
    psi_ref(i,1)=ref_eul(i,1);
end
%% Calculate Attitude using Accelerometer, Magnetometer,and Gyroscope Measurements  
    for i = 1:length(Acc)
        % First Equ.
        phi1_acc(i,1)	=	atan2(acc_y(i),sqrt(acc_x(i)^2 + acc_z(i)^2));
        theta1_acc(i,1)	=	atan2(-acc_x(i),sqrt(acc_y(i)^2 + acc_z(i)^2));
        % Second Equ.
        phi2_acc(i,1)   =   atan2(acc_y(i),acc_z(i));
        theta2_acc(i,1) =   atan2(-acc_x(i),(acc_y(i)*sin(phi2_acc(i))+acc_z(i)*cos(phi2_acc(i))));

        % Third Equ.
        theta3_acc(i,1) =   atan2(-acc_x(i),1);
        
        psi_mag(i,1)        =   deg2rad(90) + atan2(mag_z(i)*sin(theta1_acc(i,1)) - mag_y(i)*cos(theta1_acc(i,1)), mag_x(i)*cos(phi2_acc(i,1)) + sin(phi2_acc(i,1))*(mag_y(i)*sin(theta1_acc(i,1)) + mag_z(i)*cos(theta1_acc(i,1))));

    end


%% Calculate Phi, Theta, and Psi using  Measurements  
    phi_gyro    =   zeros(length(Acc),1);
    theta_gyro  =   zeros(length(Acc),1);
    psi_gyro    =   zeros(length(Acc),1);
    
    phi_CF1   = zeros(length(Acc),1);
    phi_CF2   = zeros(length(Acc),1);
    phi_CF1(1,1)   =   phi2_acc(1,1);
    phi_CF2(1,1)   =   phi2_acc(1,1);
    theta_CF1 = zeros(length(Acc),1);
    theta_CF2 = zeros(length(Acc),1);
    theta_CF3 = zeros(length(Acc),1);
    theta_CF1(1,1) =   theta1_acc(1,1);
    theta_CF2(1,1) =   theta1_acc(1,1);
    theta_CF3(1,1) =   theta1_acc(1,1);

    psi_CF    = zeros(length(Acc),1);
    psi_CF(1,1)   =   psi_mag(1,1);
    alpha = 0.3;
     for i = 2:length(Acc)
        p   =   gyro_x(i);
        q   =   gyro_y(i);
        r   =   gyro_z(i);

        phi     =   phi_CF2(i-1);
        theta    =  theta_CF1(i-1);
        psi     =   psi_mag(i-1);
        phi_gyro(i,1)   =     phi     +   dt*(p + (sin(phi)*tan(theta)*q) + (cos(phi)*tan(theta)*r));
        theta_gyro(i,1) =     theta   +   dt*((cos(phi)*q) - (sin(phi)*r));
        psi_gyro(i,1)   =     psi     +   dt*(q*(sin(phi)/cos(theta)) + r*(cos(phi)/cos(theta)));

        phi_CF1(i,1)	= ((1-alpha) * phi_gyro(i))     + (alpha * phi1_acc(i));
        phi_CF2(i,1)	= ((1-alpha) * phi_gyro(i))     + (alpha * phi2_acc(i));
        
        theta_CF1(i,1)	= ((1-alpha) * theta_gyro(i))	+ (alpha * theta1_acc(i));
        theta_CF2(i,1)	= ((1-alpha) * theta_gyro(i))	+ (alpha * theta2_acc(i));
        theta_CF3(i,1)	= ((1-alpha) * theta_gyro(i))	+ (alpha * theta3_acc(i));

        psi_CF(i,1)     = ((1-alpha) * psi_gyro(i))	+ (alpha * psi_mag(i));
    end
    
    figure(1)
    plot(t,rad2deg(phi_CF1),t,rad2deg(phi_CF2),t,rad2deg(phi_ref))
    legend('CF1','CF2','Ref')
    title('Phi CF')
    RMSE_phi_CF1 = sqrt(mean( (rad2deg(phi_CF1)-rad2deg(phi_ref)).^2 ));
    RMSE_phi_CF2 = sqrt(mean( (rad2deg(phi_CF2)-rad2deg(phi_ref)).^2 ));
    
    figure(2)
    plot(t,rad2deg(theta_CF1),t,rad2deg(theta_CF2),t,rad2deg(theta_CF3),t,rad2deg(theta_ref))
    legend('CF1','CF2','CF3','Ref')
    title('Theta CF')
    RMSE_theta_CF1 = sqrt(mean( (rad2deg(theta_CF1)-rad2deg(theta_ref)).^2 ));
    RMSE_theta_CF2 = sqrt(mean( (rad2deg(theta_CF2)-rad2deg(theta_ref)).^2 ));
    RMSE_theta_CF3 = sqrt(mean( (rad2deg(theta_CF3)-rad2deg(theta_ref)).^2 ));

    figure(3)
    plot(t,rad2deg(psi_CF),t,rad2deg(psi_ref))
    legend('CF','Ref')
    title('Psi CF')
    RMSE_psi_CF = sqrt(mean( (rad2deg(psi_CF)-rad2deg(psi_ref)).^2 ));
    
    RMSE_CF = [RMSE_phi_CF1, RMSE_phi_CF2,0
                RMSE_theta_CF1, RMSE_theta_CF2, RMSE_theta_CF3
                RMSE_psi_CF,0 ,0]
end

 

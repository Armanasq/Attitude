%% Atitude Determination From Acc and Mag
% inputs:   Accelerometer Measurements  Acc:  n x 3 matric
%           Magnetometer Measurements   Mag:  n x 3 matric
%           Reference Quaternion        Quat: n x 4 matric
%           Sample Rate Frequency       fs:   1 x 1 scaler (Hz)
function [RMSE_acc] = AccMagCal(Acc,Mag,Quat,fs)

for i = 1:length(Acc)
    acc_x(i,1) = Acc(i,1)/norm(Acc(i,:));
    acc_y(i,1) = Acc(i,2)/norm(Acc(i,:));
    acc_z(i,1) = Acc(i,3)/norm(Acc(i,:));
    mag_x(i,1) = Mag(i,1)/norm(Mag(i,:));
    mag_y(i,1) = Mag(i,2)/norm(Mag(i,:));
    mag_z(i,1) = Mag(i,3)/norm(Mag(i,:));
end
%% normalizing the measurements
[A,b,expmfs] = magcal(Mag); % calibration coefficients
%expmfs
Mag = (Mag-b)*A; % calibrated data
for i = 1:length(Acc)
    mag_x(i,1) = Mag(i,1);
    mag_y(i,1) = Mag(i,2);
    mag_z(i,1) = Mag(i,3);
end
% calculate the time
    dt      =   1/fs;
    t       =   0:dt:(length(Acc) - 1)*dt;
%% Calculate Attitude Reference
ref_eul = quat2eul(Quat);
phi_ref(:,1) = ref_eul(:,3);
theta_ref(:,1)=ref_eul(:,2);
psi_ref(:,1)=ref_eul(:,1);

%% Calculate Attitude using Accelerometer and Magnetometer Measurements  

% Accelerometer
    for i = 1:length(Acc)
        % First Equ.
        phi1_acc(i,1)	=	atan2(acc_y(i),sqrt(acc_x(i)^2 + acc_z(i)^2));
        theta1_acc(i,1)	=	atan2(-acc_x(i),sqrt(acc_y(i)^2 + acc_z(i)^2));
        % Second Equ.
        phi2_acc(i,1)   =   atan2(acc_y(i),acc_z(i));
        theta2_acc(i,1) =   atan2(-acc_x(i),(acc_y(i)*sin(phi2_acc(i))+acc_z(i)*cos(phi2_acc(i))));

        % Third Equ.
        theta3_acc(i,1) =   atan2(-acc_x(i),1);
        
        psi(i,1)        =   deg2rad(90) + atan2(mag_z(i)*sin(theta1_acc(i,1)) - mag_y(i)*cos(theta1_acc(i,1)), mag_x(i)*cos(phi2_acc(i,1)) + sin(phi2_acc(i,1))*(mag_y(i)*sin(theta1_acc(i,1)) + mag_z(i)*cos(theta1_acc(i,1))));

    end
     
    
    figure(1)
    plot(t,rad2deg(phi1_acc),t,rad2deg(phi2_acc),t,rad2deg(phi_ref),LineWidth = 1.2)
    legend('Phi1','Phi2', 'Ref')
    title('Phi Acc')
    RMSE_phi1_acc = sqrt(mean( (rad2deg(phi1_acc)-rad2deg(phi_ref)).^2 ));
    RMSE_phi2_acc = sqrt(mean( (rad2deg(phi2_acc)-rad2deg(phi_ref)).^2 ));

    figure(2)
    plot(t,rad2deg(theta1_acc),t,rad2deg(theta2_acc),t,rad2deg(theta3_acc),t,rad2deg(theta_ref),LineWidth = 1.2)
    legend('Theta1','Theta2','Theta3', 'Ref')
    title('Theta Acc')
    RMSE_theta1_acc = sqrt(mean( (rad2deg(theta1_acc)-rad2deg(theta_ref)).^2 ));
    RMSE_theta2_acc = sqrt(mean( (rad2deg(theta2_acc)-rad2deg(theta_ref)).^2 ));
    RMSE_theta3_acc = sqrt(mean( (rad2deg(theta3_acc)-rad2deg(theta_ref)).^2 ));

    figure(3)
    plot(t,rad2deg(psi),t,rad2deg(psi_ref),LineWidth = 1.2)
    legend('Psi', 'Ref')
    title('Psi Mag')
    RMSE_psi = sqrt(mean( (rad2deg(psi)-rad2deg(psi_ref)).^2 ));

    RMSE_acc = [RMSE_phi1_acc,RMSE_phi2_acc,0
            RMSE_theta1_acc, RMSE_theta2_acc, RMSE_theta3_acc
            RMSE_psi, 0 , 0]
    
end
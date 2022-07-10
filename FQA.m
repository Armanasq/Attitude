%% Factored Quaternion Algorithm 
% inputs:   Sensor Measurements         input:  n x 6 matric        [Acc(n x 3) , Mag(n x 3)]
%           Reference Quaternion        Quat:   n x 4 matric
%           height of the loc.          height  1 x 1 scaler (m)
%           location informations       latitude  1 x 1
%                                       longitude 1 x 1
%           recording date              year, month, day
%           Sample Rate Frequency       fs:   1 x 1 scaler (Hz)
function [RMSE_FQA]= FQA(acc, mag,Quat,fs, height,latitude,longitude, year,mon,day)
%% calculate the time
dt      =   1/fs;
t       =   0:dt:(length(acc) - 1)*dt;
%% Calculate Attitude Reference
ref_eul         =   quat2eul(Quat);
phi_ref(:,1)    =   ref_eul(:,3);
theta_ref(:,1)  =   ref_eul(:,2);
psi_ref(:,1)    =   ref_eul(:,1);
%% DET_R, Det_Att
[XYZ, ~, ~, I, ~] = wrldmagm(height,latitude,longitude, decyear(year,mon,day),'2015v2');
%% Referece vector data
Mx = cosd(I); %XYZ(1,1)/norm(XYZ);
My = 0;%XYZ(2,1)/norm(XYZ);
Mz = -sind(I); %XYZ(3,1)/norm(XYZ);


for i = 1:length(acc)
    acc(i,1) = acc(i,1)/norm(acc(i,:));
    acc(i,2) = acc(i,2)/norm(acc(i,:));
    acc(i,3) = acc(i,3)/norm(acc(i,:));
% Elevation Quaternion
    s_theta = acc(i,1);
    c_theta = sqrt(1 - s_theta^2);
    
    s_theta2 = sign(s_theta) * sqrt( (1-c_theta)/2 );
    c_theta2 = sqrt( (1+c_theta)/2 );

    q_elv =[c_theta2, 0, s_theta2, 0 ];
    %q_elv = q_elv/norm(q_elv);

    if c_theta == 0
        s_phi = 0;
        c_phi = 0;
    else
        s_phi = acc(i,2)/c_theta;
        c_phi = acc(i,3)/c_theta;
    end

    if c_phi >= 1
        c_phi = 1;
    elseif c_phi <= -1
        c_phi = -1;
    end
    if (sign(s_phi)==0) OR (c_phi ==-1 )
        s_phi2 =  1 * sqrt( (1 - c_phi)/2 ) ; 
    else
        s_phi2 =  sign(s_phi) * sqrt( (1 - c_phi)/2 ) ; 
    end
    c_phi2 = sqrt( (1 + c_phi)/2 );
    q_roll = ([c_phi2, s_phi2, 0 ,0]);
    q_roll = q_roll/norm(q_roll);
    
    mag(i,1) = mag(i,1)/norm(mag(i,:));
    mag(i,2) = mag(i,2)/norm(mag(i,:));
    mag(i,3) = mag(i,3)/norm(mag(i,:));
    mag_temp = [mag(i,1), mag(i,2), mag(i,3)];
    bm = [0, mag_temp];
    q_rCe = quatmultiply(quatinv(q_roll),quatinv(q_elv));
    q_bm_rCe = quatmultiply(bm,q_rCe);
    q_r_bm_rCe = quatmultiply(q_roll,q_bm_rCe);
    em = quatmultiply(q_elv, q_r_bm_rCe);
    em = em/norm(em);
    mx = em(1,2);
    my = em(1,3);
    nx = Mx/norm(Mx,My);
    ny = My/norm(Mx,My);
    psi = [mx my;-my mx]*[nx;ny];
    c_psi = psi(1,1);
    s_psi = psi(2,1);
    s_psi2 = sign(s_psi)*sqrt( (1-c_psi)/2 );
    c_psi2 = sqrt( (1+c_psi)/2 );
    q_az = [c_psi2,0,0,s_psi2];
    %q_az = q_az/norm(q_az);

    q_e_r = quatmultiply(q_elv,q_roll);
    q_e_r = q_e_r/norm(q_e_r);
    

      q_opt = quatmultiply(q_az,q_e_r);
      q_opt = q_opt/norm(q_opt);

    att_quat(i,:)   =    q_opt;
    fixTerm = [0.7071068 0 0  0.7071];
    %att_quat(i,:)   =    quatmultiply(fixTerm,q_opt); 
%     quat_err        =    quatmultiply(att_quat(i,:),quatinv(Quat(i,:)));
%     err(i,:)        =    2*acos(sqrt(quat_err(1,1)^2 + quat_err(1,4)^2));
end
% RMSE_quest          =    rad2deg(sqrt(mean(err.^2)))
% % Convert Quat_Est to Euler Angles 
est_eul = quat2eul(att_quat);

phi_est(:,1)    = est_eul(:,3);
theta_est(:,1)  = est_eul(:,2);
psi_est(:,1)    = est_eul(:,1);
% 
% figure(1)
% plot(t,rad2deg(phi_est),t,rad2deg(phi_ref))
% RMSE_phi_dav = sqrt(mean( (rad2deg(phi_est)-rad2deg(phi_ref)).^2 ));
% 
% figure(2)
% plot(t,rad2deg(theta_est), t ,rad2deg(theta_ref))
% RMSE_theta_dav = sqrt(mean( (rad2deg(theta_est)-rad2deg(theta_ref)).^2 ));
% 
% figure(3)
% plot(t,rad2deg(psi_est),t,rad2deg(psi_ref))
% RMSE_psi_dav = sqrt(mean( (rad2deg(psi_est)-rad2deg(psi_ref(1000,1))).^2 ));

figure(4)
subplot(3,1,1)
plot(psi_est+pi/2)
title('yaw')
subplot(3,1,2)
plot(-theta_est)
title('pitch')
subplot(3,1,3)
plot(phi_est)
title('roll')


figure(5)
subplot(3,1,1)
plot(psi_ref)
title('yaw')
subplot(3,1,2)
plot(theta_ref)
title('pitch')
subplot(3,1,3)
plot(phi_ref)
title('roll')
end
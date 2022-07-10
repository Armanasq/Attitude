%% Davenport’s q-Method
function [RMSE_dav]= davenport(acc, mag,Quat,fs, height,latitude,longitude, year,mon,day)
%% calculate the time
dt      =   1/fs;
t       =   0:dt:(length(acc) - 1)*dt;
%% Calculate Attitude Reference
ref_eul         =   quat2eul(Quat);
phi_ref(:,1)    =   ref_eul(:,3);
theta_ref(:,1)  =   ref_eul(:,2);
psi_ref(:,1)    =   ref_eul(:,1);
%% DET_R, Det_Att
[XYZ, H, D, I, F] = wrldmagm(height,latitude,longitude, decyear(year,mon,day),'2015v2');
%% Referece vector data
Mx = cosd(I); %XYZ(1,1)/norm(XYZ);
My = 0;%XYZ(2,1)/norm(XYZ);
Mz = -sind(I); %XYZ(3,1)/norm(XYZ);
M  = [Mx;My;Mz];

gx = 0;
gy = 0;
gz = 1;
G  = [gx;gy;gz];

for i = 1:length(acc)
    acc(i,1) = acc(i,1)/norm(acc(i,:));
    acc(i,2) = acc(i,2)/norm(acc(i,:));
    acc(i,3) = acc(i,3)/norm(acc(i,:));
    acc_temp = [acc(i,1);acc(i,2);acc(i,3)];
    
    mag(i,1) = mag(i,1)/norm(mag(i,:));
    mag(i,2) = mag(i,2)/norm(mag(i,:));
    mag(i,3) = mag(i,3)/norm(mag(i,:));
    mag_tem  = [mag(i,1);mag(i,2); mag(i,3)];
    
    b1 = acc_temp*G';
    b2 = mag_tem*M';

    B = 0.5*(b1+b2);
    S = B+transpose(B);
    z = [B(2,3)-B(3,2)
         B(3,1)-B(1,3)
         B(1,2)-B(2,1)];
    K = [trace(B) transpose(z)
            z     (S-trace(B)*eye(3))];
    
    [V,D] = eig(K,'nobalance');      
    [d,id] = max(diag(D));
    att_quat(i,:) = quatmultiply([0.7071068 0 0  0.7071],V(:,id)'); 
    quat_err = quatmultiply(att_quat(i,:),quatinv(Quat(i,:)));
    err(i,:)      = 2*acos(sqrt(quat_err(1,1)^2 + quat_err(1,4)^2));
end
RMSE_dav = rad2deg(sqrt(mean(err.^2)))
%% Convert Quat_Est to Euler Angles 
est_eul = quat2eul(att_quat);

phi_est(:,1)    = est_eul(:,3);
theta_est(:,1)  = est_eul(:,2);
psi_est(:,1)    = est_eul(:,1);

figure(1)
plot(t,rad2deg(phi_est),t,rad2deg(phi_ref))
RMSE_phi_dav = sqrt(mean( (rad2deg(phi_est)-rad2deg(phi_ref)).^2 ));

figure(2)
plot(t,rad2deg(theta_est), t ,rad2deg(theta_ref))
RMSE_theta_dav = sqrt(mean( (rad2deg(theta_est)-rad2deg(theta_ref)).^2 ));

figure(3)
plot(t,rad2deg(psi_est),t,rad2deg(psi_ref))
RMSE_psi_dav = sqrt(mean( (rad2deg(psi_est)-rad2deg(psi_ref)).^2 ));

end
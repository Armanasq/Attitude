%% Triaxial Attitude Determination (TRIAD) Method
% H. D. Black, "A passive system for determining the attitude of a satellite," AIAA journal, 1964.
% inputs:   Sensor Measurements         input:  n x 6 matric        [Acc(n x 3) , Mag(n x 3)]
%           Reference Quaternion        Quat:   n x 4 matric
%           height of the loc.          height  1 x 1 scaler (m)
%           location informations       latitude  1 x 1
%                                       longitude 1 x 1
%           recording date              year, month, day
%           Sample Rate Frequency       fs:   1 x 1 scaler (Hz)
%  TRIAD([input1(:,1:3),input1(:,7:9)],output1,fs, 44, 52.51491,13.3268,2019,7,16);
function [RMSE_TRIAD] = TRIAD(input,Quat,fs, height,latitude,longitude, year,mon,day)
%% calculate the time
dt      =   1/fs;
t       =   0:dt:(length(input) - 1)*dt;
%% Calculate Attitude Reference
ref_eul = quat2eul(Quat);
phi_ref(:,1) = ref_eul(:,3);
theta_ref(:,1)=ref_eul(:,2);
psi_ref(:,1)=ref_eul(:,1);
%% DET_R, Det_Att
[XYZ, H, D, I, F] = wrldmagm(height,latitude,longitude, decyear(year,mon,day),'2015v2');


%% Referece vector data
Mx = cosd(I); %XYZ(1,1)/norm(XYZ);
My = 0;%XYZ(2,1)/norm(XYZ);
Mz = -sind(I); %XYZ(3,1)/norm(XYZ);
gx = 0;
gy = 0;
gz = -1;
M  = [Mx;My;Mz];
G  = [gx;gy;gz];

%% Ref TRIAD
R1 = G;
R2 = cross(G,M)/norm(cross(G,M));
R3 = cross(R1,R2)/norm(cross(R1,R2));
R_ref = [R1,R2,R3];

%% Observation
ax = input(:,1);
ay = input(:,2);
az = input(:,3);

mx = input(:,4);
my = input(:,5);
mz = input(:,6);
mag = [mx,my,mz];
%% Mag Cal
for i = 1: length(input)
    x(1,i) = mx(i,1)'/norm(mag(i,:));%
    y(1,i) = my(i,1)'/norm(mag(i,:));
    z(1,i) = mz(i,1)'/norm(mag(i,:));
end
D = [x(:),y(:),z(:)];
[A,b,expmfs] = magcal(D); % calibration coefficients
%expmfs
C = (D-b)*A; % calibrated data
%C = D;
att_quat = zeros(length(input),4);
    for i = 1:length(input)
        ax(i,1) = ax(i,1)/norm(input(i,1:3));
        ay(i,1) = ay(i,1)/norm(input(i,1:3));
        az(i,1) = az(i,1)/norm(input(i,1:3));
        acc_temp = [ax(i,1);ay(i,1);az(i,1)];
        mag_temp = [C(i,1);C(i,2);C(i,3)];
        S1 = acc_temp/norm(acc_temp);
        S2 = cross(acc_temp,mag_temp)/norm(cross(acc_temp,mag_temp));
        S3 = cross(S1,S2)/norm(cross(S1,S2));
        R_Meas = [S1,S2,S3];
        Attitude = R_ref*R_Meas';
        Attitude = Attitude*[0 1 0;1 0 0; 0 0 -1];
        att_quat(i,:) = dcm2quat(Attitude)/norm(dcm2quat(Attitude));
        quat_err = quatmultiply(att_quat(i,:),quatinv(Quat(i,:)));
        err(i,:)      = 2*acos(sqrt(quat_err(1,1)^2 + quat_err(1,4)^2));
    end
    RMSE_TRIAD = rad2deg(sqrt(mean(err.^2)))

%% Convert Quat_Est to Euler Angles 
est_eul = quat2eul(att_quat);

phi_est(:,1) = est_eul(:,3);
theta_est(:,1)=est_eul(:,2);
psi_est(:,1)=est_eul(:,1);

figure(1)
plot(t,rad2deg(phi_est),t,rad2deg(phi_ref))

figure(2)
plot(t,rad2deg(theta_est), t ,rad2deg(theta_ref))

figure(3)
plot(t,rad2deg(psi_est),t,rad2deg(psi_ref))
end
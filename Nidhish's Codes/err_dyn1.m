clear; clc

T = 5;
dt = 0.001;
t = [0:dt:T]';
N = length(t);

R = zeros(3,3,N);
omg = zeros(3,N);
del = zeros(1,N);
del_dot = zeros(1,N);
euler = zeros(3,N);

% initial conditions
%R0 = eye(3);
axang0 = [0 0 1 pi/2];
R0 = axang2rotm(axang0);
omg0 = [1*pi/2;0;0*pi/2];
del0 = 0;
del_dot0 = 0;

R(:,:,1) = R0;
omg(:,1) = omg0;
del(:,1) = del0;
del_dot(:,1) = del_dot0;

R1=R0; omg1=omg0; del1=del0; del_dot1=del_dot0;
u1 = zeros(4,1);
u_all = zeros(4,N);
del_d1 = 0;
del_d = zeros(1,N);

% des rot-traj
input_freq = 2*pi*1;
input_amp = 20*pi/180;
des_ax = [1;0;1]; des_ax = des_ax/norm(des_ax);
des_ax = des_ax.*ones(3,N);
axis_ang_d = [des_ax;input_amp*sin(input_freq*t)'];
omega_d = input_amp*input_freq*cos(input_freq*t)'.*des_ax;
alpha_d = -input_amp*(input_freq)^2*sin(input_freq*t)'.*des_ax;
alpha_dd = -input_amp*(input_freq)^3*cos(input_freq*t)'.*des_ax;
alpha_ddd = input_amp*(input_freq)^4*sin(input_freq*t)'.*des_ax;
euler_d = zeros(3,N);

% model params
tau_d_max = 0.8; % (N-m)
Ixx = 1.1116548*1e-2; % wing-1,2 inertia about com (kg-m^2)
Iyy = 1.3604647*1e-2; Izz = 2.2756692*1e-2; %
I0 = diag([Ixx,Iyy,Izz]);
l = 0.42; % dist b/w wings (m)
T0 = 5.886; % thrust (N)

% control params
% k_R = diag([1.7, 1.7, 0.07]); k_omg = diag([0.279,0.279,0.02]);
% eta_n = [2;2;0.5]; zeta = [1;1;1]; % natural freq and damping
k_R = I0.*diag((2*pi*[1,1,1]).^2); k_omg = I0.*diag((2*2*pi*[1,1,1]));
eta_n = 2*pi*[1;1;1]; zeta = [1;1;1];
D = diag(2*zeta.*eta_n); K = diag(eta_n.^2);
A = [zeros(3), eye(3); -K, -D]; B = [zeros(3);eye(3)];
M = [0;0;0]; Md = M;
Mx = [M;Md];

for i=2:N
    
    del1 = 0;
    
    omg_des = omega_d(:,i);
    omg_desd = alpha_d(:,i);
    omg_desdd = alpha_dd(:,i);
    omg_desddd = alpha_ddd(:,i);
    %control
    I = 2*[Ixx,0,0;
       0, cos(del1)^2*Iyy + sin(del1)^2*Izz,0;
       0,0,cos(del1)^2*Izz + sin(del1)^2*Iyy];
    I_inv = diag(1./[I(1,1),I(2,2),I(3,3)]);
    R_d = axang2rotm(axis_ang_d(:,i)');
    Re = R_d'*R1;
    M = Mx(1:3);  Md = Mx(4:6);
    omgd = I_inv*(M-hat(omg1)*I*omg1);
    omgdd = I_inv*(Md - hat(omgd)*I*omg1 - hat(omg1)*I*omgd);
    e_omg = omg1 - Re'*omg_des;
    e_omgd = omgd - Re'*omg_desd + hat(e_omg)*Re'*omg_des;
    e_omgdd = omgdd + hat(e_omg)*Re'*omg_desd - Re'*omg_desdd ...
        + hat(e_omgd)*Re'*omg_des - hat(e_omg)^2*Re'*omg_des + hat(e_omg)*Re'*omg_desd;
    e_R = vee(Re-Re')/2;
    e_Rd = vee(Re*hat(e_omg) + hat(e_omg)*Re')/2;
    e_Rdd = vee(Re*hat(e_omg)^2 + hat(e_omgd)*Re' - hat(e_omg)^2*Re')/2;
    M_des = -k_R*e_R - k_omg*e_omg + Re'*(I*omg_desd + hat(omg_des)*I*omg_des);
    M_desd = -k_R*e_Rd - k_omg*e_omgd - hat(e_omg)*Re'*(I*omg_desd + hat(omg_des)*I*omg_des)...
        + Re'*(I*omg_desdd + hat(omg_desd)*I*omg_des + hat(omg_des)*I*omg_desd);
    M_desdd = -k_R*e_Rdd - k_omg*e_omgdd + (hat(e_omg)^2*Re'-hat(e_omgd)*Re')*(I*omg_desd + hat(omg_des)*I*omg_des)...
        -hat(e_omg)*Re'*(I*omg_desdd + hat(omg_desd)*I*omg_des + hat(omg_des)*I*omg_desd)...
        + Re'*(I*omg_desddd + hat(omg_desdd)*I*omg_des + 2*hat(omg_desd)*I*omg_desd + hat(omg_des)*I*omg_desdd);
    u = M_desdd + D*M_desd + K*M_des;
    Mx = Mx + (A*Mx + B*u)*dt;
%     Mz = Mx(3); Mzd = Mx(6);
%     v1 = -(u(3) + 2*Mz*Mzd/(l*T0) - D(3,3)*Mzd - K(3,3)*Mz)*l*T0/(Mz^2 + (l*T0)^2);
%     del_des = atan(-M_des(3)/(l*T0));
%     del_desd = -l*T0*M_desd(3)/(M_des(3)^2 + (l*T0)^2);
%     del_desdd = (-l*T0*M_desdd(3) - 2*M_des(3)*M_desd(3))/(M_des(3)^2 + (l*T0)^2);
%     v1 = del_desdd - 2*2*pi*2*(del_dot1 - del_desd) - (2*pi*2)^2*(del1 - del_des);
%     c1 = sin(2*del1)*(Iyy-Izz)*(omg1(2)^2 - omg1(3)^2)/(2*Ixx);
%     u1(2) = Ixx*(c1+v1);
%     u1(1) = Mx(1);
%     u1(3) = Mx(2)/l;
%     del_d1 = Mx(3);
%     del_d1 = (abs(del_d1)>pi/6)*sign(del_d1)*pi/6 + (abs(del_d1)<=pi/6)*del_d1;
%     u1(2) = -4*sin(del1 - del_d1) - 0.4*del_dot1;
%     u1(2) = (abs(u1(2))>tau_d_max)*sign(u1(2))*tau_d_max + (abs(u1(2))<=tau_d_max)*u1(2);
    u1(1:3) = Mx(1:3);
    u1(4) = 5.886; %Thrust (N)
    u_all(:,i) = u1;
%     del_d(i) = del_d1;
    
    %dynamics
    %[R1,del1,omg1,del_dot1] = sw_biplane_model(R1,del1,omg1,del_dot1,u1,dt);
    omg1 = omg1 + dt*I_inv*(Mx(1:3) - hat(omg1)*I*omg1);
    R1 = R1*expm(hat(omg1)*dt);
    
    R(:,:,i) = R1;
    del(i) = del1;
    omg(:,i) = omg1;
    del_dot(:,i) = del_dot1;
    euler(:,i) = rotm2eul(R1);
    euler_d(:,i) = rotm2eul(R_d);
    
end

%% plot

figure
plot(t,omg); grid on;
title('\omega');

% figure
% plot(t,del); hold on
% plot(t,del_dot); grid on
% legend('\delta','\delta_{dot}')
% title('\delta')
% 
% figure
% plot(t,del_d*180/pi); grid on; hold on
% plot(t,del*180/pi); grid on; hold on
% legend('\delta_{des}','\delta')
% title('\delta_{des}')

figure
plot(t,u_all(1:3,:)); grid on
title('u')

figure
plot(t,euler*180/pi); grid on; hold on
plot(t,euler_d*180/pi);
legend('yaw','pit','roll')
title('Euler Angle')

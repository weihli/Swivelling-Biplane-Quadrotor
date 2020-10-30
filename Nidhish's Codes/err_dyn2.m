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
euler_d = zeros(3,N);

% model params
tau_d_max = 0.8; % (N-m)
Ixx = 1.1116548*1e-2; % wing-1,2 inertia about com (kg-m^2)
Iyy = 1.3604647*1e-2; Izz = 2.2756692*1e-2; %
I0 = diag([Ixx,Iyy,Izz]);
l = 0.42; % dist b/w wings (m)
T0 = 5.886; % thrust (N)

% control params
k_R = I0.*diag((2*pi*[1,1,1]).^2); k_omg = I0.*diag((2*2*pi*[1,1,1]));

for i=2:N
    
    del1 = 0;
    
    omg_des = omega_d(:,i);
    omg_desd = alpha_d(:,i);
    
    %control
    I = 2*[Ixx,0,0;
       0, cos(del1)^2*Iyy + sin(del1)^2*Izz,0;
       0,0,cos(del1)^2*Izz + sin(del1)^2*Iyy];
    I_inv = diag(1./[I(1,1),I(2,2),I(3,3)]);
    R_d = axang2rotm(axis_ang_d(:,i)');
    Re = R_d'*R1;
    e_omg = omg1 - Re'*omg_des;
    e_R = vee(Re-Re')/2;
    M_des = -k_R*e_R - k_omg*e_omg + Re'*(I*omg_desd + hat(omg_des)*I*omg_des);
%     M_des = -k_R*e_R - k_omg*e_omg + hat(omg1)*I*omg1 ...
         - I*(hat(omg1)*Re'*omg_des - Re'*omg_desd);
    u1(1:3) = M_des;
    u1(4) = 5.886; %Thrust (N)
    u_all(:,i) = u1;
    
    %dynamics
    omg1 = omg1 + dt*I_inv*(M_des - hat(omg1)*I*omg1);
    R1 = R1*expm(hat(omg1)*dt);
    
    R(:,:,i) = R1;
    del(i) = del1;
    omg(:,i) = omg1;
    del_dot(:,i) = del_dot1;
    euler(:,i) = rotm2eul(R1)';
    euler_d(:,i) = rotm2eul(R_d)';
    
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
subplot(3,1,1);plot(t,euler(1,:)*180/pi); grid on; hold on
plot(t,euler_d(1,:)*180/pi); legend('yaw', 'yaw_d')
subplot(3,1,2);plot(t,euler(2,:)*180/pi); grid on; hold on
plot(t,euler_d(2,:)*180/pi); legend('pit', 'pit_d')
subplot(3,1,3);plot(t,euler(3,:)*180/pi); grid on; hold on
plot(t,euler_d(3,:)*180/pi); legend('roll', 'roll_d')
title('Euler Angle')

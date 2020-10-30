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
del(:,1) = del0;
omg(:,1) = omg0;
del_dot(:,1) = del_dot0;

R1=R0; del1=del0; omg1=omg0; del_dot1=del_dot0;
u1 = zeros(4,1);
u = zeros(4,N);
del_d1 = 0;
del_d = zeros(1,N);

% des rot-traj
input_freq = 2*pi*1;
input_amp = 20*pi/180;
des_ax = [1*ones(1,N);0*ones(1,N);0*ones(1,N)];
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

for i=2:N
    
    %control
    R_d = axang2rotm(axis_ang_d(:,i)');
    Re = R_d'*R1;
    e_R = vee((Re-Re'))/2;
    M_d = -1.7*e_R - 0.2794*omg1;
    u1(1) = M_d(1);
    u1(3) = M_d(2)/0.42;
    del_d1 = 0.07*e_R(3) + 0.02*omg1(3);
    del_d1 = (abs(del_d1)>pi/6)*sign(del_d1)*pi/6 + (abs(del_d1)<=pi/6)*del_d1;
    u1(2) = -4*sin(del1 - del_d1) - 0.4*del_dot1;
    u1(2) = (abs(u1(2))>tau_d_max)*sign(u1(2))*tau_d_max + (abs(u1(2))<=tau_d_max)*u1(2);
    u1(4) = 5.886; %Thrust (N)
    u(:,i) = u1;
    del_d(i) = del_d1;
    
    %dynamics
    [R1,omg1,del1,del_dot1] = sw_biplane_model(R1,omg1,del1,del_dot1,u1,dt);
    
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

figure
plot(t,del); hold on
plot(t,del_dot); grid on
legend('\delta','\delta_{dot}')
title('\delta')

figure
plot(t,del_d*180/pi); grid on; hold on
plot(t,del*180/pi); grid on; hold on
legend('\delta_{des}','\delta')
title('\delta_{des}')

figure
plot(t,u(1:3,:)); grid on
title('u')

figure
plot(t,euler*180/pi); grid on; hold on
plot(t,euler_d*180/pi);
legend('yaw','pit','roll')
title('Euler Angle')

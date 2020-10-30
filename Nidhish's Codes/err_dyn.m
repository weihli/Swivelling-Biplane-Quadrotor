clear; clc

T = 5;
dt = 0.001;
t = [0:dt:T]';
N = length(t);

R = zeros(3,3,N);
omg = zeros(3,N);
Me = zeros(3,N);
euler = zeros(3,N);

axang0 = [0 0 1 pi/2];
R0 = axang2rotm(axang0);
omg0 = [1*pi/2;0;0*pi/2];
M0 = [0;0;0];

R(:,:,1) = R0;
omg(:,1) = omg0;
Me(:,1) = M0;

Ixx = 1.1116548*1e-2; % wing-1,2 inertia about com (kg-m^2)
Iyy = 1.3604647*1e-2; Izz = 2.2756692*1e-2; %
I = diag([Ixx,Iyy,Izz]);
I_inv = diag([1/Ixx,1/Iyy,1/Izz]);

k_R = I.*diag((2*pi*[1,1,1]).^2); k_omg = I.*diag((2*2*pi*[1,1,1]));
eta_n = 2*pi*[1;1;1]; zeta = [1;1;1]; % natural freq and damping
D = diag(2*zeta.*eta_n); K = diag(eta_n.^2);
A = [zeros(3), eye(3); -K, -D]; B = [zeros(3);eye(3)];
M = [1;0;0]; Md = M;
Mx = [M;Md];

omg1 = omg0;
R1 = R0;


for i=2:N
    e_R = vee(R1-R1')/2;
    e_omg = omg1;
    
    Mx = Mx + A*Mx*dt;
    M = Mx(1:3);
    omg1 = omg1 + I_inv*(-k_R*e_R - k_omg*e_omg + M)*dt;
    R1 = R1*expm(hat(omg1)*dt);
    
    Me(:,i) = M;
    omg(:,i) = omg1;
    euler(:,i) = rotm2eul(R1);
end

figure
plot(t,euler*180/pi);
legend('yaw','pitch','roll')
title('Euler')

figure
plot(t,omg);
title('\omega')

figure
plot(t,Me);
title('M')
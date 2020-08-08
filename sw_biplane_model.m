function  [R_1,omg_1,del_1,del_dot_1] = sw_biplane_model(R,omg,del,del_dot,u,h)

R_del = [1,0,0;
        0,cos(del),-sin(del);
        0,sin(del),cos(del)];

% Rotational drag
KD = 0*diag([0.001,0.001,0.025]);
    
% inertia matrix about center of mass
Ixx = 1.1116548*1e-2; Iyy = 1.3604647*1e-2; Izz = 2.2756692*1e-2;
I = 2*[Ixx,0,0;
       0, cos(del)^2*Iyy + sin(del)^2*Izz,0;
       0,0,cos(del)^2*Izz + sin(del)^2*Iyy];
I_inv = diag([1/I(1,1),1/I(2,2),1/I(3,3)]);
I_dot = [0,0,0;0,1,0;0,0,-1]*sin(2*del)*(Izz-Iyy)*del_dot;

l = 0.42; %m

tau_m = u(1);
tau_d = u(2);
T_d = u(3);
T_m = u(4);

M = 2*[tau_m;l*T_d*cos(del);-l*T_m*sin(del)] - norm(omg)*KD*omg;

%omg_d = I_inv*(-I_dot*omg - hat(omg)*I*omg + M);
% omg_1 = omg + h*omg_d;
k1 = h*I_inv*(-I_dot*omg - hat(omg)*I*omg + M);
k2 = h*I_inv*(-I_dot*(omg+k1/2) - hat(omg+k1/2)*I*(omg+k1/2) + M);
k3 = h*I_inv*(-I_dot*(omg+k2/2) - hat(omg+k2/2)*I*(omg+k2/2) + M);
k4 = h*I_inv*(-I_dot*(omg+k3) - hat(omg+k3)*I*(omg+k3) + M);
omg_1 = omg + (k1+2*k2+2*k3+k4)/6;

del_ddot = (2*tau_d - sin(2*del)*(Iyy-Izz)*(omg(2)^2 - omg(3)^2))/(2*Ixx);
del_dot_1 = del_dot + h*del_ddot;
omg_del_1 = [del_dot_1;0;0];

R_1 = R*expm(hat(omg_1)*h);
R_del_1 = R_del*expm(hat(omg_del_1)*h);
del_1 = atan2(R_del_1(3,2),R_del_1(2,2));

end
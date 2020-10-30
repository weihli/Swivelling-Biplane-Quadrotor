clear; clc

T = 5;
dt = 0.001;
t = [0:dt:T]';
N = length(t);

% states
R = zeros(3,3,N);
omg = zeros(3,N);
del = zeros(1,N);
del_dot = zeros(1,N);
euler = zeros(3,N);

% initial conditions
axang0 = [0 0 1 pi/2];
R0 = axang2rotm(axang0);
omg0 = [0;0;0*pi/2];
del0 = 0;
del_dot0 = 0;

R(:,:,1) = R0;
omg(:,1) = omg0;
del(:,1) = del0;
del_dot(:,1) = del_dot0;

R1=R0; omg1=omg0; del1=del0; del_dot1=del_dot0;
u1 = zeros(4,1);
u_all = zeros(4,N);
F_all = zeros(4,N);
del_d1 = 0;
del_d = zeros(1,N);

% des rot-traj
input_freq = 2*pi*1;
input_amp = 40*pi/180;
des_ax = [0;1;1]; des_ax = des_ax/norm(des_ax);
des_ax = des_ax.*ones(3,N);
axis_ang_d = [des_ax;input_amp*sin(input_freq*t)'];
omega_d = input_amp*input_freq*cos(input_freq*t)'.*des_ax;
alpha_d = -input_amp*(input_freq)^2*sin(input_freq*t)'.*des_ax;
alpha_dd = -input_amp*(input_freq)^3*cos(input_freq*t)'.*des_ax;
alpha_ddd = input_amp*(input_freq)^4*sin(input_freq*t)'.*des_ax;
euler_d = zeros(3,N);

% model params
tau_d_max = 0.8; % (N-m)
% Ixx = 0.95*1.1116548*1e-2; % wing-1,2 inertia about com (kg-m^2)
% Iyy = 0.95*1.3604647*1e-2; Izz = 1.05*2.2756692*1e-2; %
Ixx = 1.1116548*1e-2; % wing-1,2 inertia about com (kg-m^2)
Iyy = 1.3604647*1e-2; Izz = 2.2756692*1e-2; %
I0 = diag([Ixx,Iyy,Izz]);
l = 0.42; % dist b/w wings (m)
L = 0.6; % dist b/w motors on same wing (m)
T0 = 5.886; % thrust (N)
A = [-L/2, L/2, 0,0; 0,0,-L/2, L/2; 1,1,0,0; 0,0,1,1];
B = [-1,1,0,0;1,1,0,0;0,0,1,1;0,0,-1,1];
G = inv(A)*B; G_inv = inv(G);
f_min = 0.2145; f_min = ones(4,1)*f_min; % min motor force (N)
f_max = 6.74; f_max = ones(4,1)*f_max;  % max motor force (N)

% control params
k_R = I0.*diag((2*pi*[4,4,1.16]).^2); k_omg = I0.*diag((2*2*pi*[1.5,1.5,1]));
eta_n = 2*pi*[6;6;4]; zeta = [1;1;1];
eta_del = 2*pi*4; zeta_del = 1;
D = diag(2*zeta.*eta_n); K = diag(eta_n.^2);
A = [zeros(3), eye(3); -K, -D]; B = [zeros(3);eye(3)];
M = zeros(3,1); Md = M;
Mx = [M;Md];
Mz = 0; Mzd = 0;
u1d = zeros(4,1);

for i=2:N
    
    omg_des = omega_d(:,i);
    omg_desd = alpha_d(:,i);
    omg_desdd = alpha_dd(:,i);
    omg_desddd = alpha_ddd(:,i);
    %control
    I = diag(2*[Ixx,Iyy,Izz]);
    I_inv = diag(1./[I(1,1),I(2,2),I(3,3)]);
    R_d = axang2rotm(axis_ang_d(:,i)');
    Re = R_d'*R1;
    M = Mx(1:3);  Md = Mx(4:6);
    M(3) = Mz; Md(3) = Mzd;
    omgd = I_inv*M;
    omgdd = I_inv*Md;
    e_omg = omg1 - Re'*omg_des;
    e_omgd = omgd - Re'*omg_desd;
    e_omgdd = omgdd - Re'*omg_desdd;
    e_R = vee(Re-Re')/2;
    e_Rd = vee(Re*hat(e_omg) + hat(e_omg)*Re')/2;
    e_Rdd = vee(Re*hat(e_omgd) + hat(e_omgd)*Re')/2;
%     M_des = -k_R*e_R - k_omg*e_omg + Re'*I*omg_desd;
%     M_desd = -k_R*e_Rd - k_omg*e_omgd + Re'*I*omg_desdd;
%     M_desdd = -k_R*e_Rdd - k_omg*e_omgdd + Re'*I*omg_desddd;
    M_des = -k_R*e_R - k_omg*e_omg + Re'*(I*omg_desd + hat(omg_des)*I*omg_des);
    M_desd = -0*k_R*e_Rd - 0*k_omg*e_omgd...
        + Re'*(I*omg_desdd + hat(omg_desd)*I*omg_des + hat(omg_des)*I*omg_desd);
    M_desdd = -0*k_R*e_Rdd - 0*k_omg*e_omgdd...
        + Re'*(I*omg_desddd + hat(omg_desdd)*I*omg_des + 2*hat(omg_desd)*I*omg_desd + hat(omg_des)*I*omg_desdd);
    u = M_desdd + D*M_desd + K*M_des;
    Mx = Mx + (A*Mx + B*u)*dt;
%     Mz = Mx(3); Mzd = Mx(6);
    Mz = -2*l*T0*tan(del1); Mzd = -2*l*T0*sec(del1)^2*del_dot1;
%     v1 = -(u(3) + 2*Mz*Mzd/(l*T0) - D(3,3)*Mzd - K(3,3)*Mz)*l*T0/(Mz^2 + (l*T0)^2);
    del_des = atan2(-M_des(3),2*l*T0);
    del_desd = -2*l*T0*M_desd(3)/(M_des(3)^2 + (2*l*T0)^2);
    del_desdd = -2*l*T0*M_desdd(3)/(M_des(3)^2 + (2*l*T0)^2) ...
        + 4*l*T0*M_des(3)*M_desd(3)^2/(M_des(3)^2 + (2*l*T0)^2)^2;
    v1 = del_desdd - 2*eta_del*(del_dot1 - del_desd) - eta_del^2*(del1 - del_des);
    c1 = sin(2*del1)*(Iyy-Izz)*(omg1(2)^2 - omg1(3)^2)/(2*Ixx);
    u1(2) = Ixx*(0*c1+v1);  % tau_delta
%     u1(1) = Mx(1)/2;  % tau_m
%     u1(3) = Mx(2)/(2*l*cos(del1));  % T_delta
    u1(1) = M_des(1)/2;  % tau_m
    u1(3) = M_des(2)/(2*l);%*cos(del1));  % T_delta
%     del_d1 = Mx(3);
%     del_d1 = (abs(del_d1)>pi/6)*sign(del_d1)*pi/6 + (abs(del_d1)<=pi/6)*del_d1;
    u1(2) = (abs(u1(2))>tau_d_max)*sign(u1(2))*tau_d_max + (abs(u1(2))<=tau_d_max)*u1(2);
%     u1(4) = 5.886; %Thrust (N)
    u1(4) = T0;%/cos(del1); % T_m
    f1 = G*u1;
    f1 = (f1<f_min).*f_min + (f1>f_min).*f1;
    f1 = (f1>f_max).*f_max + (f1<f_max).*f1;
    u1 = G_inv*f1;
    F_all(:,i) = f1;
    u_all(:,i) = u1;
    % actuator dynamics
    u1d = u1d + dt/0.015*(u1-u1d);
%     u1d = u1;
    %dynamics
    [R1,omg0,del1,del_dot1] = sw_biplane_model(R1,omg0,del1,del_dot1,u1d,dt);
    % sensor noise
    omg1 = omg0 + 0*0.075*rand(3,1);
    
    R(:,:,i) = R1;
    omg(:,i) = omg1;
    del(i) = del1;
    del_dot(:,i) = del_dot1;
    del_d(i) = del_des;
    euler(:,i) = rotm312eul(R1);
    euler_d(:,i) = rotm312eul(R_d);
    
end

euler(:,1) = rotm312eul(R0);
euler_d(:,1) = rotm312eul(axang2rotm(axis_ang_d(:,1)'));

%% plot

font_size = 9;
lfont_size = 9;
fig_dim = [0 0 14 13.5];
figure('Units','centimeters','PaperPosition',fig_dim,'PaperPositionMode','manual',...
    'PaperSize',[14 13.5])
set(gcf,'DefaultLineLineWidth',1)
subplot(3,1,1)
plot(t,omg(1,:)); grid on; hold on
plot(t,omega_d(1,:),'-.')
legend('\omega_x','\omega_{dx}')
set(gca,'FontSize',font_size,'FontName','Times','FontWeight','normal')
legend('boxoff'); legend('Orientation','horizontal')
legend('Location','best')
set(gca,'XTickLabel',[])
title('\omega');
subplot(3,1,2)
plot(t,omg(2,:)); grid on; hold on
plot(t,omega_d(2,:),'-.')
legend('\omega_y','\omega_{dy}')
set(gca,'FontSize',font_size,'FontName','Times','FontWeight','normal')
legend('boxoff'); legend('Orientation','horizontal')
legend('Location','best')
set(gca,'XTickLabel',[])
subplot(3,1,3)
plot(t,omg(3,:)); grid on; hold on
plot(t,omega_d(3,:),'-.')
legend('\omega_z','\omega_{dz}')
set(gca,'FontSize',font_size,'FontName','Times','FontWeight','normal')
legend('boxoff'); legend('Orientation','horizontal')
legend('Location','best')
xlabel('Time (s)')
%savefig('fig/omg')
%print -dpdf fig/omg.pdf

figure
plot(t,u_all(1:3,:)); grid on
title('u')

fig_dim = [0 -1 14 22.5];
figure('Units','centimeters','PaperPosition',fig_dim,'PaperPositionMode','manual',...
    'PaperSize',[14 21.5])
set(gcf,'DefaultLineLineWidth',1)
subplot(5,1,1)
plot(t,euler(1,:)*180/pi); grid on; hold on
plot(t,euler_d(1,:)*180/pi,'-.');
set(gca,'FontSize',font_size,'FontName','Times','FontWeight','normal')
ylabel('\psi (deg)')
legend('\psi','\psi_d')
legend('boxoff'); legend('Orientation','horizontal')
legend('Location','best')
set(gca,'XTickLabel',[])
title('312 Euler Angles')

subplot(5,1,2)
plot(t,euler(2,:)*180/pi); grid on; hold on
plot(t,euler_d(2,:)*180/pi,'-.');
set(gca,'FontSize',font_size,'FontName','Times','FontWeight','normal')
ylabel('\phi (deg)')
legend('\phi','\phi_d')
legend('boxoff'); legend('Orientation','horizontal')
legend('Location','best')
set(gca,'XTickLabel',[])

subplot(5,1,3)
plot(t,euler(3,:)*180/pi); grid on; hold on
plot(t,euler_d(3,:)*180/pi,'-.');
set(gca,'FontSize',font_size,'FontName','Times','FontWeight','normal')
ylabel('\theta (deg)')
legend('\theta','\theta_d')
legend('boxoff'); legend('Orientation','horizontal')
legend('Location','best')
set(gca,'XTickLabel',[])

subplot(5,1,4)
plot(t,del*180/pi); grid on; hold on
plot(t,del_d*180/pi,'-.');
set(gca,'FontSize',font_size,'FontName','Times','FontWeight','normal')
ylabel('\delta (deg)')
legend('\delta','\delta_d')
legend('boxoff'); legend('Orientation','horizontal')
legend('Location','best')
title('Swivel Angle')

subplot(5,1,5)
plot(t,F_all(1,:),t,F_all(2,:),':k',t,F_all(3,:),'-.r',t,F_all(4,:),'--m'); grid on
set(gca,'FontSize',font_size,'FontName','Times','FontWeight','normal')
ylabel('Force (N)')
legend('f_1','f_2','f_3','f_4')
legend('boxoff'); legend('Orientation','horizontal')
legend('Location','best')
title('Motor Force')

xlabel('Time (s)')
% savefig('fig/angle')
% print -dpdf fig/angle.pdf

% fig_dim = [0 0 14 4.5];
% figure('Units','centimeters','PaperPosition',fig_dim,'PaperPositionMode','manual',...
%     'PaperSize',[14 4.5])
% set(gcf,'DefaultLineLineWidth',1)
% plot(t,F_all(1,:),t,F_all(2,:),':k',t,F_all(3,:),'-.r',t,F_all(4,:),'--m'); grid on
% set(gca,'FontSize',font_size,'FontName','Times','FontWeight','normal')
% legend('f_1','f_2','f_3','f_4')
% legend('boxoff'); legend('Orientation','horizontal')
% legend('Location','best')
% xlabel('Time (s)')
% title('Motor Force (N)');
% savefig('fig/motor')
% print -dpdf fig/motor.pdf
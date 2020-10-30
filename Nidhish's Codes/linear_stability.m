Ixx = 1.1116548*1e-2; % wing-1,2 inertia about com (kg-m^2)
Iyy = 1.3604647*1e-2; Izz = 2.2756692*1e-2; %
I0 = diag([Ixx,Iyy,Izz]);

I_inv = diag([1/I0(1,1),1/I0(2,2),1/I0(3,3)]);

fi = 4; % inner loop freq (Hz)
fo = 1; % outer loop freq (Hz)

k_R = I0.*diag((2*pi*fo*[1,1,1]).^2); k_omg = I0.*diag((2*2*pi*fo*[1,1,1]));
eta_n = 2*pi*fi*[1;1;1]; zeta = [1;1;1];
D = diag(2*zeta.*eta_n); K = diag(eta_n.^2);

I = eye(3); Z = zeros(3);

P = diag([1,2,3]); B = zeros(3);
%R = eye(3);
[vec1, val1] = eig(P);
R = expm(pi*hat(vec1(:,3)));
for i = 1:3
B = B - 0.5*hat(I(:,i))*P*R*hat(I(:,i));
end

S = [Z, I, Z, Z;
    -I_inv*k_R*B, -I_inv*k_omg, I_inv, Z;
    Z, Z, Z, I;
    Z, Z, -K, -D];

[vec, val] = eig(S);
diag(val)
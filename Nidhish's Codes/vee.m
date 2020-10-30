function omega = vee(omega_hat)
omega = zeros(3,1);
omega(1,1) = omega_hat(3,2);
omega(2,1) = omega_hat(1,3);
omega(3,1) = omega_hat(2,1);
end
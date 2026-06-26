function PV = PV_inverse(x,eps0)
% principal value for a 1/x
% This is usually done for the real part of 1/(omega_nm + i*eps0), and the
% imaginary part is actually \pi * delta_funct(omega_nm, eps0)
PV = x ./ (x.^2+eps0^2);
end
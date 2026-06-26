function morb = calc_orb_BerryConn(nw,vmn,omn,ksi)
% This is based on Xie, X.C. et al. arXiv:2108.03928
morb = zeros(nw,3);
A1 = 1i*vmn(:,:,1) .* omn ./ (omn.^2 + ksi^2);
A2 = 1i*vmn(:,:,2) .* omn ./ (omn.^2 + ksi^2);
A3 = 1i*vmn(:,:,3) .* omn ./ (omn.^2 + ksi^2);
v1 = vmn(:,:,1);
v2 = vmn(:,:,2);
v3 = vmn(:,:,3);
m12 = 0.5*real(v1*A2 - A1*v2);
m23 = 0.5*real(v2*A3 - A2*v3);
m31 = 0.5*real(v3*A1 - A3*v1);
morb(:,3) = diag(m12)/2/pi;
morb(:,2) = diag(m31)/2/pi;
morb(:,1) = diag(m23)/2/pi;
end
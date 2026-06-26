function morbit ...
    = calc_orb_ModernTheory_v2(energy, vmn, omn, nw, EF, ksi)
morbit = zeros(nw,3);
omplusn = energy*ones(1,nw) + ones(nw,1)*energy';
mu = EF * eye(nw,nw);
fac1 = 2*mu - omplusn;
vx = vmn(:,:,1);
vy = vmn(:,:,2);
vz = vmn(:,:,3);
rx = -1i*vx .* omn ./ (omn.^2 + ksi^2);
ry = -1i*vy .* omn ./ (omn.^2 + ksi^2);
rz = -1i*vz .* omn ./ (omn.^2 + ksi^2);

% These are in unit of e^2*V*angstrom^2/hbar
morbit(:,1) = diag(imag((ry.*fac1) * rz - rz * (ry.*fac1)));
morbit(:,2) = diag(imag((rz.*fac1) * rx - rx * (rz.*fac1)));
morbit(:,3) = diag(imag((rx.*fac1) * ry - ry * (rx.*fac1)));

morbit = morbit * 0.13124; % To convert e/hbar * eV * [r^2] to mu_B (e*hbar/2/m_e)

end
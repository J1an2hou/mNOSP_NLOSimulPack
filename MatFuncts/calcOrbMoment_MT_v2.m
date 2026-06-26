function morbit ...
    = calcOrbMoment_MT_v2(energy, vmn, omn, nw, EF, ksi)
morbit = zeros(nw,3);
omplusn = energy*ones(1,nw) + ones(nw,1)*energy';
mu = EF * eye(nw,nw);
fac1 = (2*mu - omplusn);
vx = vmn(:,:,1);
vy = vmn(:,:,2);
rx = -1i*vx .* omn ./ (omn.^2 + ksi^2);
ry = -1i*vy .* omn ./ (omn.^2 + ksi^2);
% fac2z = -imag(rx.*ry.' - ry.*rx.');
% morbit(:,3) = diag(fac1 * fac2z);
morbit(:,3) = -diag(imag((rx.*fac1) * ry - ry * (rx.*fac1)));

end
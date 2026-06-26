function morbit ...
    = calc_orb_ModernTheory(Hk, energy, wavefunct, vmn, omn, nw, EF, ksi)
dudk = gradu(wavefunct, vmn, omn, nw, ksi);
morbit = zeros(nw,3);
egvalues = diag(energy);
mu = eye(nw) * EF;
m12 = 2 * dudk(:,:,1)' * (Hk + egvalues - 2*mu) * dudk(:,:,2);
m23 = 2 * dudk(:,:,2)' * (Hk + egvalues - 2*mu) * dudk(:,:,3);
m31 = 2 * dudk(:,:,3)' * (Hk + egvalues - 2*mu) * dudk(:,:,1);
morbit(:,1) = imag(diag(m23));
morbit(:,2) = imag(diag(m31));
morbit(:,3) = imag(diag(m12)); % in e/(2\hbar)

end
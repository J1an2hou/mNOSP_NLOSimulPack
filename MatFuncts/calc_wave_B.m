function wave_B = calc_wave_B(wavefunct, nw, morb, omn, Bfield, ksi)
% morb_mn = < m | L | n >, denoted as Lmn
% wavefunct_n = wavefunct_n + B * \sum_m (Lmn/omn) * wavefunct_m
% See Eq. (5) in https://arxiv.org/pdf/2402.05241.pdf
% Note that this is only fine with small B, as we don't make wavefunctions
% perpendicular, just normalized

mu_b2eV = 5.78e-5; % from Tesla * mu_B to eV
L_o = zeros(nw,nw,3);
L_o(:,:,1) = morb(:,:,1) .* omn ./ (omn.^2+ksi^2);
L_o(:,:,2) = morb(:,:,2) .* omn ./ (omn.^2+ksi^2);
L_o(:,:,3) = morb(:,:,3) .* omn ./ (omn.^2+ksi^2);
B1 = Bfield(1)*mu_b2eV;
B2 = Bfield(2)*mu_b2eV;
B3 = Bfield(3)*mu_b2eV;
wB = wavefunct + B1*wavefunct * L_o(:,:,1) + B2*wavefunct * L_o(:,:,3) ...
    + B3*wavefunct * L_o(:,:,3);
wave_B = wB ./ sqrt(sum(wB .* conj(wB)));


end

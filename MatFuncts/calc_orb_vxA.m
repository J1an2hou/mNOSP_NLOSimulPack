function morb_interband = calc_orb_vxA(nw,vmn,omn,ksi)
% This is based on Yang, Shengyuan et al. PRL 132, 056301 (2024)
% Interband < m | L | n > = e/hbar/2*((v_ml x A_ln) + (v_nn x A_mn)) 
morb_interband = zeros(nw,nw,3);
A1 = -1i*vmn(:,:,1) .* omn ./ (omn.^2 + ksi^2);
A2 = -1i*vmn(:,:,2) .* omn ./ (omn.^2 + ksi^2);
A3 = -1i*vmn(:,:,3) .* omn ./ (omn.^2 + ksi^2);
v1 = vmn(:,:,1);
v2 = vmn(:,:,2);
v3 = vmn(:,:,3);
diagv1 = real(diag(vmn(:,:,1)));
diagv2 = real(diag(vmn(:,:,2)));
diagv3 = real(diag(vmn(:,:,3)));
vnn1 = ones(nw,1) * diagv1';
vnn2 = ones(nw,1) * diagv2';
vnn3 = ones(nw,1) * diagv3';

m1 = (v2 * A3 - v3 * A2) + (vnn2 .* A3 - vnn3 .* A2);
m2 = (v3 * A1 - v1 * A3) + (vnn3 .* A1 - vnn1 .* A3);
m3 = (v1 * A2 - v2 * A1) + (vnn1 .* A2 - vnn2 .* A1);
% Alternatively, one can use
% vmnplus1 = ones(nw,1) * diagv1' + diagv1 * ones(1,nw);
% vmnplus2 = ones(nw,1) * diagv2' + diagv2 * ones(1,nw);
% vmnplus3 = ones(nw,1) * diagv3' + diagv3 * ones(1,nw);
% m1 = 1i*((A2.*omn) * A3 - (A3.*omn) * A2) + (vmnplus2 .* A3 - vmnplus3 .* A2);
% m2 = 1i*((A3.*omn) * A1 - (A1.*omn) * A3) + (vmnplus3 .* A1 - vmnplus1 .* A3);
% m3 = 1i*((A1.*omn) * A2 - (A2.*omn) * A1) + (vmnplus1 .* A2 - vmnplus2 .* A1);
% End of the alternative algorithm

% These are in unit of e^2*V*angstrom^2/hbar/2

morb_interband(:,:,1) = m1;
morb_interband(:,:,2) = m2;
morb_interband(:,:,3) = m3;

%Convert to mu_B
morb_interband = morb_interband * 0.13124;

end
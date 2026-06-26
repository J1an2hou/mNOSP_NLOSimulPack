function m = calc_orb_vxA_v2(nw,vmn,rmn,dir)
% This is based on Yang, Shengyuan et al. PRL 132, 056301 (2024)
% Interband < m | L | n > = e/hbar/2*((v_ml x A_ln) + (v_nn x A_mn)) 
A1 = rmn(:,:,1);
A2 = rmn(:,:,2);
A3 = rmn(:,:,3);
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

% These are in unit of e^2*V*angstrom^2/hbar/2, then convert into mu_B
m1 = m1 * 0.13124;
m2 = m2 * 0.13124;
m3 = m3 * 0.13124;

if nargin == 3
    m = zeros(nw, nw, 3);
    m(:,:,1) = m1;
    m(:,:,2) = m2;
    m(:,:,3) = m3;

% Until now this is the usual case
elseif nargin == 4
    if dir == 1
        m = m1;
    elseif dir == 2
        m = m2;
    elseif dir == 3
        m = m3;
    else
        error('Check your spin direction, which is not 1 - 3')
    end
else
    error('Check your input arguments for calc_orb_vxA_v2 function!')
end

end
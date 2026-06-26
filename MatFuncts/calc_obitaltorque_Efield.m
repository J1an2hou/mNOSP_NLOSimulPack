function [tau_cp, tau_on_band] = calc_obitaltorque_Efield(vmn, rmn, energy, ChemPot, nw, L_dir, E_dir, ksi)
% L_dir and E_dir are directions for orbital and electric field
% For example, for 2D materials, L_dir should be 3, and E_dir can be 1 or 2
% It takes from Eq. (4) in PRB 110, 035427 (2024), but it seems that the
% second term (three band process) is absolutely zero under \epsilon_ijk
ncp = length(ChemPot);
tau_on_band = zeros(nw, ncp);
tau_cp = zeros(1, ncp);
[a, b] = inv_Voigt1D(L_dir);
ImQb = imag(rmn(:,:,b) .* (rmn(:,:,E_dir).'));
ImQa = imag(rmn(:,:,a) .* (rmn(:,:,E_dir).'));
va = diag(real(vmn(:,:,a)));
vb = diag(real(vmn(:,:,b)));
va_sum = va*ones(1,nw) + ones(nw,1)*va';
vb_sum = vb*ones(1,nw) + ones(nw,1)*vb';
fac1 = diag(ImQb * va_sum - ImQa * vb_sum); % in (nw, 1) dimension
% r_Edir = rmn(:,:,E_dir);
% vb_OD = vmn(:,:,b) - diag(diag(vmn(:,:,b)));
% va_OD = vmn(:,:,a) - diag(diag(vmn(:,:,a)));
% vba_anticom = vmn(:,:,b) * va_OD + vmn(:,:,a) * vb_OD;
% vab_anticom = vmn(:,:,a) * vb_OD + vmn(:,:,b) * va_OD;
% fac2 = diag(real((r_Edir .* omn ./ (omn.^2 + ksi^2)) ...
%     * (vba_anticom-vab_anticom)));
for icp = 1:ncp
    mu = energy - ChemPot(icp);
    fnn = 1./(1+exp(mu/ksi));
    tau_on_band(:, icp) = fac1 ; %- fac2;
    tau_cp(icp) = sum(fnn .* tau_on_band(:, icp));
end

end
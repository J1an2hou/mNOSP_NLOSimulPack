function [dcdt_sp,chi_FS_sp,chi_dis_sp,dcdt_nl,chi_FS_nl,chi_dis_nl] ...
    = LinearEdel_sp_Neel(v,o,en,sp1,sp2,sp3,nl1,nl2,nl3,nw,cp,EF,direct,ksi)
% We evaluate all three kinds of E-M coupling
% The first one is dominant in metallic phase, Fermi surface contributed
% It is evaluated by d chi / dt, so we didn't multiply relaxation lifetime
% The second is the intrinsic semiconductor term, contributed by Fermi sea
% The last one scales with disorder energy ksi, hence can be viewed as
% disorder contributed term. In most cases, this should be small.
e_charge = -1;
ncp = length(cp);
ndir = size(direct,1);
dcdt_sp = zeros(ndir, ncp);
chi_FS_sp = zeros(ndir, ncp);
chi_dis_sp = zeros(ndir, ncp);
dcdt_nl = zeros(ndir, ncp);
chi_FS_nl = zeros(ndir, ncp);
chi_dis_nl = zeros(ndir, ncp);

for icp = 1:ncp
mu = EF + cp(icp);
f = 1./(1+exp((en-mu)/ksi));
% dfde = ksi ./ ((energy-mu).^2+ksi^2)/pi;
fmn = f*ones(1,nw) - ones(nw,1)*f';
deltaFS = delta_funct(en - mu, ksi);
facFS = fmn .* (o.^2 - ksi^2) ./ ((o.^2 + ksi^2).^2);
facdis = fmn .* o ./ ((o.^2 + ksi^2).^2);
for idir = 1:ndir
    a = direct(idir,1);
    b = direct(idir,2);
    va = v(:,:,a);
    if b == 1
        sp = sp1; nl = nl1;
    elseif b == 2
        sp = sp2; nl = nl2;
    elseif b == 3
        sp = sp3; nl = nl3;
    else
        error('Wrong spin or Neel vector direction')
    end
    dcdt_sp(idir, icp) = sum(diag(real(sp)) .* diag(real(va)) .* deltaFS);
    dcdt_nl(idir, icp) = sum(diag(real(nl)) .* diag(real(va)) .* deltaFS);
    chi_FS_sp(idir, icp) = imag(trace((facFS.*sp)*va - (facFS.*va)*sp)/2);
    chi_FS_nl(idir, icp) = imag(trace((facFS.*nl)*va - (facFS.*va)*nl)/2);
    chi_dis_sp(idir, icp) = -ksi*real(trace((facdis.*sp)*va + (facdis.*va)*sp));
    chi_dis_nl(idir, icp) = -ksi*real(trace((facdis.*nl)*va + (facdis.*va)*nl));
end

end
% converting units, 241.8 is converting h to ps (241.8 = 1/0.004135)
dcdt_sp = dcdt_sp * e_charge / 10 * 241.8; % in unit of mu_B * nm / V /ps
dcdt_nl = dcdt_nl * e_charge / 10 * 241.8; % in unit of mu_B * nm / V /ps
chi_FS_sp = chi_FS_sp * e_charge / 10; % in unit of mu_B * nm / V
chi_FS_nl = chi_FS_nl * e_charge / 10; % in unit of mu_B * nm / V
chi_dis_sp = chi_dis_sp * e_charge / 10; % in unit of mu_B * nm / V
chi_dis_nl = chi_dis_nl * e_charge / 10; % in unit of mu_B * nm / V
end
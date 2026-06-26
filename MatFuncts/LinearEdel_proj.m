function [dchidt,chi_FS,chi_dis] ...
    = LinearEdel_proj(v,o,en,p,sp1,sp2,sp3,nw,cp,EF,direct,ksi)
% We evaluate all three kinds of linear Edelstein effect
% The first one is dominant in metallic phase, Fermi surface contributed
% It is evaluated by d chi / dt, so we didn't multiply relaxation lifetime
% This is a time-reversal even contribution
% The second is the intrinsic semiconductor term, contributed by Fermi sea
% This is a time-reversal odd contribution
% The last one scales with disorder energy ksi, hence can be viewed as
% disorder contributed term. In most cases, this should be small (also
% time-reversal even contribution)

e_charge = -1;
ncp = length(cp);
ndir = size(direct,1);
nwl = size(p,3);
dchidt = zeros(ndir, ncp, nwl);
chi_FS = zeros(ndir, ncp, nwl);
chi_dis = zeros(ndir, ncp, nwl);

for icp = 1:ncp
mu = EF + cp(icp);
f = 1./(1+exp((en-mu)/ksi));
% dfde = ksi ./ ((energy-mu).^2+ksi^2)/pi;
fmn = f*ones(1,nw) - ones(nw,1)*f';
deltaFS = delta_funct(en - mu, ksi);
facFS = fmn .* (o.^2 - ksi^2) ./ (o.^2 + ksi^2).^2;
facdis = fmn .* o ./ (o.^2 + ksi^2).^2;
for idir = 1:ndir
    a = direct(idir,1);
    b = direct(idir,2);
    va = v(:,:,a);
    if b == 1
        sp = sp1;
    elseif b == 2
        sp = sp2;
    elseif b == 3
        sp = sp3;
    else
        error('Wrong spin or Neel vector direction')
    end
    for iwl = 1:nwl
    pj = p(:,:,iwl);
    sp_proj = comm(sp, pj, 1)/2;
    dchidt(idir, icp, iwl) = sum(diag(real(sp_proj)) .* diag(real(va)) .* deltaFS);
    chi_FS(idir, icp, iwl) = imag(trace((facFS.*sp_proj)*va - (facFS.*va)*sp_proj)/2);
    chi_dis(idir, icp, iwl) = -ksi*real(trace((facdis.*sp_proj)*va + (facdis.*va)*sp_proj));
    end
end
end
% converting units, 241.8 is converting h to ps (241.8 = 1/0.004135)
dchidt = dchidt * e_charge / 10 * 241.8; % in unit of mu_B * nm / V /ps
chi_FS = chi_FS * e_charge / 10; % in unit of mu_B * nm / V
chi_dis = chi_dis * e_charge / 10; % in unit of mu_B * nm / V
end
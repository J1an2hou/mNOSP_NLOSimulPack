function [tork_odd,tork_even] ...
    = SOTtorkance_proj(v,o,en,p,tor,nw,cp,EF,direct,ksi)
% We evaluate SOT effect - electric field induced torque
% The torkance factor can be found in https://arxiv.org/pdf/1305.4873
% It contains generally two parts
% T-odd part is from fermi surface contribution
% At the clean limit, it is tork ~ < n | T | n > < n | v | n > *FermiSurf
% The T even part is akin to Berry curvature
% At the clean limit, tork ~ Im ( Tnm * vmn / (En-Em)^2 )
% Here we adopt a general formalism to evaluate them

% The unit is e * angstrom

ncp = length(cp);
ndir = size(direct,1);
nwl = size(p,3);
tork_odd = zeros(ndir, ncp, nwl);
tork_even = zeros(ndir, ncp, nwl);
for icp = 1:ncp
mu = EF + cp(icp);
f = 1./(1+exp((en-mu)/ksi));
fmn = f*ones(1,nw) - ones(nw,1)*f';
deltaFS = delta_funct(en - mu, ksi);
% facFS = fmn .* (o.^2 - ksi^2) ./ (o.^2 + ksi^2).^2;
% double_FS = (deltaFS*ones(1,nw)) .* (ones(nw,1)*deltaFS')*pi^2;
% Wmn = Twoenergy_broad(en,nw,mu,ksi);

for idir = 1:ndir
    a = direct(idir,1);
    b = direct(idir,2);
    va = v(:,:,a);
    tb = tor(:,:,b);
    for iwl = 1:nwl
    pj = p(:,:,iwl);
    tor_proj = comm(tb, pj, 1)/2;
    % tv = tor_proj .* va.';
    % tork_odd(idir, icp, iwl) = trace(double_FS * real(tv))/pi;
    % tork_even(idir, icp, iwl) = trace(Wmn * imag(tv))/pi;
    tork_odd(idir,icp,iwl) = ...
        sum(real(deltaFS .* diag(tor_proj) .* diag(va)))/ksi; % This is the clean limit
    tork_even(idir,icp,iwl) = -trace(imag((fmn .* (tor_proj .* o./(o.^2+ksi^2))) ...
        * (va .* o ./ (o.^2+ksi^2)))); % This is the clean limit
    % tork_even(idir,icp,iwl) ...
    % = imag(trace((facFS.*tor_proj)*va - (facFS.*va)*tor_proj)/2); % This is another clean limit
    end
end
end

end

function Wmn_FS = Twoenergy_broad(en,nw,mu,ksi)
% We give a W function to expand the two energy bands at Fermi level
% The W = [atan((E_m - EF)/ksi)-atan((E_n-EF)/ksi)]/(E_n - E_m)^2 
%        - ksi/(E_m - E_n) * [A_n(EF)*B_m(EF) + B_n(EF)*A_m(EF)]
% One can see this equation in Eq. (S32) of PRB 111, L140409 (2025)
o = en*ones(1,nw)-ones(nw,1)*en'; % energy difference
fac1 = atan((en-mu)/ksi) * ones(1,nw) - ones(nw,1)*atan((en-mu)/ksi)';
fac1 = fac1 .* o.^2 ./ (o.^4 + ksi^4);
An = ksi ./ ((en - mu).^2 + ksi^2);
Bn = (en-mu) ./ ((en - mu).^2 + ksi^2);
fac2 = -ksi .* o ./ (o.^2 + ksi^2);
fac2 = fac2 .* (An .* Bn' + Bn .* An');
Wmn_FS = fac1 + fac2;

end
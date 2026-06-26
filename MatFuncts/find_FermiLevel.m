function [mu, Eg] = find_FermiLevel(E,Nk,ne_per_cell,ksi,Ntrial)
% The eigenenergy E is in dimension of (nkpts, nband)
% ksi is smearing factor
% ne_per_cell is the number of valence electrons, nocc in usual cases
% The output is a number that gives the uniform Fermi level energy

if nargin < 5
    Ntrial = 100; % number of trial times
end
if size(E,1) == Nk
    mu_vb = max(max(E(:,ne_per_cell)));
    mu_cb = min(min(E(:,ne_per_cell+1)));
    Eg = mu_cb - mu_vb;
elseif size(E,2) == Nk
    mu_vb = max(max(E(ne_per_cell,:)));
    mu_cb = min(min(E(ne_per_cell+1,:)));
    Eg = mu_cb - mu_vb;
else
    warning('The total kpoint number in eigenvalues does not consistent')
    warning('I will decide it by my own, but no responsibility questioned')
    Evec = E(:);
    mu_vb = min(Evec); mu_cb = max(Evec);
    Eg = 20141027;
end

if mu_cb > mu_vb
    flagrev = true;
else
    flagrev = false;
end

for itmu = 1:Ntrial
    mu_mid = 0.5*(mu_vb+mu_cb);
    ftrial = 1./(1 + exp((E - mu_mid)/ksi));
    Ntot = sum(ftrial,'all')/Nk;
    % disp(Ntot)
    % disp(['Low and High ',num2str(mu_lo),' and ', num2str(mu_hi)])
    if Ntot > ne_per_cell
        if flagrev
            mu_cb = mu_mid;
        else
            mu_vb = mu_mid;
        end
    else
        if flagrev
            mu_vb = mu_mid;
        else
            mu_cb = mu_mid;
        end
    end
    if abs(mu_cb - mu_vb) < 1e-6
        break
    end
end
mu = 0.5*(mu_vb+mu_cb);
if itmu == Ntrial
    warning('The Fermi level is not found to be precisely. Usual with care')
end

end
function H = Floquet_Recursive_SelfEnergy(R,HR,wtR,wc,nw,kpt,A0,omega,nq,alatt,blatt)
% We apply a recursive method to evaluate the self-energy of Floquet band
% engineering. In this case, one can apply an arbitrary frequency incident
% photon into the system, not limited by the high frequency limit, as
% widely used for Magnus expansion. In the first order truncation scheme,
% the self-energy can be written as
% Sigma(E) = H_{-1} / (E +1i*eta - Hk - omega) * H_{1}, and then the eigenenergy E
% is solved by self-consistent algorithm of
% (H0 + Sigma(E)) * phi = E * phi

ksi = 0.01; % a small imaginary number
if nq ~= 1
    warning('Current recursive self energy is only for first order truncation')
end
nsc = 500; % total cycle for self consistent field
thres_energy = 1e-6; % energy convergence criterion

nR = size(R,1);
FloPhase0 = FloPhase_ave(A0,R,wc,alatt,nw,nR,0);
Hk0 = buildHk(R,HR.*FloPhase0,wtR,wc,nw,nR,kpt,alatt,blatt);

FloPhaseq = FloPhase_ave(A0,R,wc,alatt,nw,nR,1);
FloPhase_q = FloPhase_ave(A0,R,wc,alatt,nw,nR,-1);
Hk1 = buildHk(R,HR.*FloPhaseq,wtR,wc,nw,nR,kpt,alatt,blatt);
Hk_1 = buildHk(R,HR.*FloPhase_q,wtR,wc,nw,nR,kpt,alatt,blatt);

E_old = eig(Hk0);
E_old = sort(real(E_old));
% self consistent evaluation
Sigma0 = zeros(nw, nw);
for isc = 1:nsc
    Sigma1 = Hk_1 / (diag(E_old) + 1i*ksi - Hk0 - Sigma0) * Hk1;
    Sigma_1 = Hk1 / (diag(E_old) + 1i*ksi - Hk0 - Sigma0) * Hk_1;
    Sigma0 = Hk1 / (diag(E_old) + 1i*ksi - Hk0 - eye(nw,nw)*omega - Sigma1) * Hk_1 ...
        +Hk_1 / (diag(E_old) + 1i*ksi - Hk0 + eye(nw,nw)*omega - Sigma_1) * Hk1;
    H = Hk0 + Sigma0;
    E_new = eig(H);
    E_new = sort(real(E_new));
    dE = norm(E_new - E_old);
    if dE < thres_energy
        break
    else
        E_old = E_new;
    end
end
if isc == nsc
    warning(['Self consistent solution for self energy not reached convergence with energy diff of ',num2str(dE)])
end

end
function Hk = TB_buildHam_onsite(tsoc,U_orbit,U_atom,m,nat,norb)
% This function is k-indpendent
% ======== Parameters ============
% tsoc: parameter
% Delta_onsite: on-site energy, not correlation
% m: magnetic exchange in Zeeman field
% ===== end Parameters ===========
nw = nat*norb*2;
Hk = zeros(nw,nw);

Hsoc11 = zeros(norb,norb);
Hsoc12 = zeros(norb,norb);
Hsoc21 = zeros(norb,norb);
Hsoc22 = zeros(norb,norb);

% on-site SOC (dxy, dx2-y2), make changes for the orbitals you selected
Hsoc11(1,2) = 2*1i* tsoc;
Hsoc11(2,1) = -2*1i* tsoc;
Hsoc22(1,2) = -2*1i* tsoc;
Hsoc22(2,1) = 2*1i* tsoc;
% Hamiltonian in basis |site-1,orb-1,up>, |site-1,orb-2,up>, ...|site-1,orb-n,up>
% |site-2,orb-1,up>, ...|site-2,orb-n,up>, ...|site-n,orb-n,up>, |site-1,orb-1,dn> ...
for iat = 1:nat
    iup1 = (iat-1)*norb+1;
    iup2 = iup1 + norb - 1;
    idn1 = iup1 + nat * norb;
    idn2 = idn1 + norb - 1;
    
    % ---------------------Zeeman field-------------------
    mx = m(iat,1)*eye(norb);
    my = m(iat,2)*eye(norb);
    mz = m(iat,3)*eye(norb);
    Hk(iup1:iup2, idn1:idn2) = Hk(iup1:iup2, idn1:idn2) + mx;
    Hk(idn1:idn2, iup1:iup2) = Hk(idn1:idn2, iup1:iup2) + mx;
    Hk(iup1:iup2, idn1:idn2) = Hk(iup1:iup2, idn1:idn2) - 1i*my;
    Hk(idn1:idn2, iup1:iup2) = Hk(idn1:idn2, iup1:iup2) + 1i*my;
    Hk(iup1:iup2, iup1:iup2) = Hk(iup1:iup2, iup1:iup2) + mz;
    Hk(idn1:idn2, idn1:idn2) = Hk(idn1:idn2, idn1:idn2) - mz;
    % ----------------End of Zeeman field----------------

    % ------------------Intrisic SOC-----------------
    Hk(iup1:iup2, iup1:iup2) = Hk(iup1:iup2, iup1:iup2) + Hsoc11;
    Hk(iup1:iup2, idn1:idn2) = Hk(iup1:iup2, idn1:idn2) + Hsoc12;
    Hk(idn1:idn2, iup1:iup2) = Hk(idn1:idn2, iup1:iup2) + Hsoc21;
    Hk(idn1:idn2, idn1:idn2) = Hk(idn1:idn2, idn1:idn2) + Hsoc22;
    % ---------------End of intrinsic SOC------------

end
% -----------------On orbital potential ----------------
Hub = diag(kron([1,1],kron(ones(1,nat),U_orbit)));
Hk = Hk + Hub;
% --------------End of on orbital potential ------------
% -----------------On atom potential ----------------
Hk = Hk +kron(eye(2,2), diag(kron(U_atom,ones(1,norb))));
% --------------End of on atom potential ------------
rdn = -6;
Hk = roundn(Hk,rdn);

% Check if they are all conjugate
if ishermitian(Hk) == 0
    error('Hk is not Hermitian')
end

end
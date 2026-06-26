function [Hk,dHdk,d2Hdk2] ...
    = TB_buildHk_hop(Vdds,Vddp,Vddd,kpt,blatt,delta1,jindex1,nc1,nat,norb)
% ========== Parameters =============
% V: Slater-Koster hopping parameters
% kpt: direct k point position
% nat: total number of sites
% nc1: first nearest neighbor number
% ========== end Parameters =========

nw = nat*norb*2;
bkpt = kpt * blatt;

% Building Hamiltonian
Hk = zeros(nw,nw);
dHdk = zeros(nw,nw,3);
d2Hdk2 = zeros(nw,nw,6);

% Hamiltonian in basis |site-1,orb-1,up>, |site-1,orb-2,up>, ...|site-1,orb-n,up>
% |site-2,orb-1,up>, ...|site-2,orb-n,up>, ...|site-n,orb-n,up>, |site-1,orb-1,dn> ...
for iat = 1:nat
    iup1 = (iat-1)*norb+1;
    iup2 = iup1 + norb - 1;
    idn1 = iup1 + nat * norb;
    idn2 = idn1 + norb - 1;
    for ic = 1:nc1
        delta_ij = squeeze(delta1(iat,ic,:));
        jat = jindex1(iat,ic);
        if jat == 0 || norm(delta_ij) < 1e-3
            continue
        end
        jup1 = (jat-1)*norb+1;
        jup2 = jup1 + norb - 1;
        jdn1 = jup1 + nat * norb;
        jdn2 = jdn1 + norb - 1;
        eikr = exp(1i*bkpt*delta_ij);
        Ht = TB_SK(Vdds,Vddp,Vddd,delta_ij,norb);
        % ---------------Hopping integral---------------
        Hk(iup1:iup2, jup1:jup2) = Hk(iup1:iup2, jup1:jup2) + Ht*eikr;
        Hk(idn1:idn2, jdn1:jdn2) = Hk(idn1:idn2, jdn1:jdn2) + Ht*eikr;
        for idim = 1:3
        dHdk(iup1:iup2, jup1:jup2, idim) ...
            = dHdk(iup1:iup2, jup1:jup2, idim) + Ht*1i*delta_ij(idim)*eikr;
        dHdk(idn1:idn2, jdn1:jdn2, idim) ...
            = dHdk(idn1:idn2, jdn1:jdn2, idim) + Ht*1i*delta_ij(idim)*eikr;
        end
        for idim = 1:6
            [a,b] = voigt(idim);
            d2Hdk2(iup1:iup2, jup1:jup2, idim) ...
                = d2Hdk2(iup1:iup2, jup1:jup2, idim) - Ht*delta_ij(a)*delta_ij(b)*eikr;
            d2Hdk2(idn1:idn2, jdn1:jdn2, idim) ...
                = d2Hdk2(idn1:idn2, jdn1:jdn2, idim) - Ht*delta_ij(a)*delta_ij(b)*eikr;
        end
        % -------------End of hopping integral---------

        % -----------Rashba Haldane SOC in 2D--------
        % Hk(iup1:iup2,jdn1:jdn2) = Hk(iup1:iup2,jdn1:jdn2) + 1i*(delta(2)+1i*delta(1))*tsoc*eikr;
        % Hk(idn1:idn2,jup1:jup2) = Hk(idn1:idn2,jup1:jup2) + 1i*(delta(2)-1i*delta(1))*tsoc*eikr;
        % for idim = 1:3
        % dHdk(iup1:iup2,jdn1:jdn2,idim) = dHdk(iup1:iup2,jdn1:jdn2,idim) ...
        %     +1i*(delta(2)+1i*delta(1))*tsoc*1i*delta(idim)*eikr;
        % dHdk(idn1:idn2,jup1:jup2,idim) = dHdk(idn1:idn2,jup1:jup2,idim) ...
        %     +1i*(delta(2)-1i*delta(1))*tsoc*1i*delta(idim)*eikr;
        % end
        % for idim = 1:6
        %     [a,b] = voigt1to2(idim);
        %     d2Hdk2(iup1:iup2,jdn1:jdn2,idim) = d2Hdk2(iup1:iup2,jdn1:jdn2,idim) ...
        %         -1i*(delta(2)+1i*delta(1))*tsoc*delta(a)*delta(b)*eikr;
        %     d2Hdk2(idn1:idn2,jup1:jup2,idim) = d2Hdk2(idn1:idn2,jup1:jup2,idim) ...
        %         -1i*(delta(2)-1i*delta(1))*tsoc*delta(a)*delta(b)*eikr;
        % end
        % --------End of Rashba Haldane SOC in 2D------
        
    end

end

rdn = -6;
Hk = roundn(Hk,rdn);
dHdk = roundn(dHdk,rdn);
d2Hdk2 = roundn(d2Hdk2,rdn);
% Check if they are all conjugate
if ishermitian(Hk) == 0
    error('Hk is not Hermitian')
end

end


function Ht = TB_SK(Vdds,Vddp,Vddd,R,norb)
% R is vector that connects site i and j
% Ht is hopping (norb x norb) integral between site i and j
% Hsoc is (norb*2 x norb*2) soc Hamiltonian
% V-prefix variables are hopping integral (ppp - Vpppi, pps - Vppsigma ...)
% U is (1 x norb) on-site energy
% tsoc is SOC strength
Ht = zeros(norb,norb);
dis = sqrt(sum(R.*R));
l = R(1)/dis;
m = R(2)/dis;
n = R(3)/dis;
Ht(1,1) = 3*l^2*m^2*Vdds+(l^2+m^2-4*l^2*m^2)*Vddp+(n^2+l^2*m^2)*Vddd; % xy - xy
Ht(1,2) = 1.5*l*m*(l^2-m^2)*Vdds+2*l*m*(m^2-l^2)*Vddp+l*m*(l^2-m^2)*Vddd/2; % xy - x2-y2
Ht(2,1) = 1.5*l*m*(l^2-m^2)*Vdds+2*l*m*(m^2-l^2)*Vddp+l*m*(l^2-m^2)*Vddd/2; % x2-y2 - xy
Ht(2,2) = 0.75*(l^2-m^2)^2*Vdds+(l^2+m^2-(l^2-m^2)^2)*Vddp+(n^2+(l^2-m^2)^2/4)*Vddd; % x2-y2 - x2-y2

end

function [Hk,dHdk,d2Hdk2] = TB_Hk_soc_chirality(tsoc,kpt,blatt,delta2,jindex,nc,nu,nat,norb)
% ==========================================================
% Intrinsic SOC , in form of i * soc * \sum_<<i,j>> nu_ij * c_i^dagger * sz * c_j
% Here, nu_ij = sign(d_ik x d_kj)_z, which is +1 for left turn, -1 for right turn
% ==========================================================
nw = nat*norb*2;
bkpt = kpt * blatt;

% Building Hamiltonian
Hk = zeros(nw,nw);
dHdk = zeros(nw,nw,3);
d2Hdk2 = zeros(nw,nw,6);
for iat = 1:nat
    iup1 = (iat-1)*norb + 1;
    iup2 = iup1 + norb - 1;
    idn1 = iup1 + nat*norb;
    idn2 = idn1 + norb - 1;

    for ic2 = 1:nc
        jat = jindex(iat,ic2);
        if jat == 0
            continue
        end
        delta = squeeze(delta2(iat,ic2,:));
        eikr  = exp(1i * (bkpt * delta));
        jup1 = (jat-1)*norb + 1;
        jup2 = jup1 + norb - 1;
        jdn1 = jup1 + nat*norb;
        jdn2 = jdn1 + norb - 1;
        nu = nu(iat,ic2);
        % Spin-up block
        Hk(iup1:iup2, jup1:jup2) = ...
            Hk(iup1:iup2, jup1:jup2) + 1i * tsoc * nu * eikr;

        % Spin-down block
        Hk(idn1:idn2, jdn1:jdn2) = ...
            Hk(idn1:idn2, jdn1:jdn2) - 1i * tsoc * nu * eikr;

        % k-derivatives
        for idim = 1:3
            dHdk(iup1:iup2, jup1:jup2, idim) = ...
                dHdk(iup1:iup2, jup1:jup2, idim) ...
                - tsoc*nu*delta(idim)*eikr;

            dHdk(idn1:idn2, jdn1:jdn2, idim) = ...
                dHdk(idn1:idn2, jdn1:jdn2, idim) ...
                + tsoc*nu*delta(idim)*eikr;
        end

        for idim = 1:6
            [a,b] = voigt(idim);
            d2Hdk2(iup1:iup2, jup1:jup2, idim) ...
                = d2Hdk2(iup1:iup2, jup1:jup2, idim) - 1i*tsoc*nu*delta(a)*delta(b)*eikr;
            d2Hdk2(idn1:idn2, jdn1:jdn2, idim) ...
                = d2Hdk2(idn1:idn2, jdn1:jdn2, idim) + 1i*tsoc*nu*delta(a)*delta(b)*eikr;
        end
    end
end

end
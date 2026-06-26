function [sppump_even,sppump_odd, orbpump_even, orbpump_odd] = SpinOrbPumpingTT(Hk,Hxc,o,p,cp,EF,orb_seq,Enk,nw,wf,direct,ksi)
% This code evaluates the spin and orbital pumping according to
% Eq. (4-5) in PRB 111, L140409 (2025)
% We breakdown the on-site contributions to the pumping from individual
% sites, and the local orbital moment is adopted.
% In this expression, we apply the Torque vs torque algorithm, so that we
% denote it as TT in the function name.

% Here spin pumping is defined as magnetization dynamics induced spin
% polarization generation, not related to spin current (which will be
% another function)

ncp = length(cp);
ndir = size(direct,1);
nwl = size(p, 3);
sppump_even = zeros(ncp, ndir, nwl, nwl);
sppump_odd = zeros(ncp, ndir, nwl, nwl);
orbpump_even = zeros(ncp, ndir, nwl, nwl);
orbpump_odd = zeros(ncp, ndir, nwl, nwl);
[lx_mat,ly_mat,lz_mat] = buildJ(orb_seq); % local orbital moment

for icp = 1:ncp
    mu = cp(icp) + EF;
    fn = 1./(1+exp((Enk-mu)/ksi));
    dfde = delta_funct(Enk-mu,ksi);
    f = fn*ones(1,nw) - ones(nw,1)*fn';
    for idir = 1:ndir
    sp_dir = direct(idir,1);
    xc_dir = direct(idir,2);
    sp_mat = Pauli_mat(nw,sp_dir);
    if sp_dir == 1
        l_mat = lx_mat;
    elseif sp_dir == 2
        l_mat = ly_mat;
    elseif sp_dir == 3
        l_mat = lz_mat;
    else
        error('Direction of angular momentum wrong')
    end
    xc_tor = Pauli_mat(nw,xc_dir);
    sptor = comm(sp_mat, Hk, -1)/1i;
    orbtor = comm(l_mat, Hk, -1)/1i;
    xctor = comm(xc_tor, Hxc, -1)/1i;
    for iwl = 1:nwl
    for jwl = 1:nwl
        pA = p(:,:,iwl);
        pB = p(:,:,jwl);
        storA = comm(sptor, pA, +1)/2;
        storAmn = wf' * storA * wf;
        ltorA = comm(orbtor, pA, +1)/2;
        ltorAmn = wf' * ltorA * wf;
        xctorB = comm(xctor, pB, +1)/2;
        xctorBmn = wf' * xctorB * wf;
        fac_sp = storAmn .* xctorBmn.';
        fac_orb = ltorAmn .* xctorBmn.';
        sppump_even(icp, idir, iwl, jwl) = pi*(dfde' * real(fac_sp) * dfde);
        sppump_odd(icp, idir, iwl, jwl) ...
            = trace((f.*o./(o.^2+ksi^2)) * (imag(fac_sp).*o./(o.^2+ksi^2)));
        orbpump_even(icp, idir, iwl, jwl) = pi*(dfde' * real(fac_orb) * dfde);
        orbpump_odd(icp, idir, iwl, jwl) ...
            = trace((f.*o./(o.^2+ksi^2)) * (imag(fac_orb).*o./(o.^2+ksi^2)));
    end
    end
    end
end

end

function s = Pauli_mat(nw, dir)
A = eye(nw/2);
sigmax = [0 1; 1 0];
sigmay = [0 -1i; 1i 0];
sigmaz = [1 0; 0 -1];
sigma0 = [1 0; 0 1];
% The basis function is |1,up>, |2,up>, ..., |1,dn>, |2,dn>, ...
s0 = kron(sigma0, A);
sx = kron(sigmax, A);
sy = kron(sigmay, A);
sz = kron(sigmaz, A);
if dir == 0
    s = s0;
elseif dir == 1
    s = sx;
elseif dir == 2
    s = sy;
elseif dir == 3
    s = sz;
else
    error('Error building Pauli matrix')
end

end
function [s2sc_even,s2oc_even,s2sc_odd,s2oc_odd] = torq2SLcurr(H,r,v,o,p,cp,EF,Enk,nw,wf,direct,ksi)
% The output of this function is spin torque induced spin current (s2sc)
% and spin torque induced orbital current (s2oc)
% They both are in size of (ncp, ndir), as we have assigned the responses
% only on the layer-indexed 3
% The inputs are:
% H is the Hamiltonian that is H_tot = H_ex + H_so
% r, v, o are the usual rmn, vmn, and omm - onn
% p is the projection matrix that for AFMs now (we may extend it later)
% this p is simply the projector, not multiplied with wavefunct (pmn)

% We give T-even and T-odd contributions separately
% In this way, it evaluates the current change over time t.
% The expressions are:
%
% s2sc_even ~ f_nm * Re(j_nm * T_mn)/(E_n-E_m) + j_nn * T_nn * df/dE
% s2sc_odd ~ f_nm * Im(j_nm * T_mn) * delta(E_n - E_m) * pi
%
% Note if electric current is considered above, then the T-odd and T-even
% swaps
ncp = length(cp);
ndir = size(direct,1);
% nwl = size(p, 3);
s2sc_even = zeros(ncp, ndir);
s2oc_even = zeros(ncp, ndir);
s2sc_odd = zeros(ncp, ndir);
s2oc_odd = zeros(ncp, ndir);
proj_A = p(:,:,1); % magnetic sublattice A
proj_B = p(:,:,2); % magnetic sublattice B
proj_C = p(:,:,3); % interested layer with spin/orbital current
proj_Cmn = wf' * proj_C * wf;

for icp = 1:ncp
    mu = cp(icp) + EF;
    fn = 1./(1+exp((Enk-mu)/ksi));
    dfde = -delta_funct(Enk-mu,ksi);
    f = fn*ones(1,nw) - ones(nw,1)*fn';
    for idir = 1:ndir
        curr_dir = direct(idir,1);
        am_dir = direct(idir,2);
        tor_dir = direct(idir,3);
        vamn = v(:,:,curr_dir);
        stor = build_Pauli(nw, tor_dir);
        stor = comm(comm(stor, H, -1), proj_A+proj_B, 0.5); % total damping from A and B
        stor_mn = wf' * stor * wf;
        sp_mn = calc_sp(nw,wf,am_dir); % tensor wvf' * S * wvf
        orb_mn = calc_orb_vxA_v2(nw, v, r, am_dir); % tensor wvf' * L * wvf
        jsp = comm(comm(vamn, sp_mn, 0.5), proj_Cmn, 0.5);
        jorb = comm(comm(vamn, orb_mn, 0.5), proj_Cmn, 0.5);
        jsp_diag = diag(real(jsp));
        jorb_diag = diag(real(jorb));
        stor_diag = diag(real(stor_mn));
        s2sc_even(icp,idir) = trace((f.*o./(o.^2+ksi^2)) * real(jsp .* stor_mn.')) ...
            + sum(jsp_diag .* stor_diag .* dfde);
        s2oc_even(icp,idir) = trace((f.*o./(o.^2+ksi^2)) * real(jorb .* stor_mn.')) ...
            + sum(jorb_diag .* stor_diag .* dfde);
        s2sc_odd(icp,idir) = -pi*trace((f.*delta_funct(o,ksi)) * imag(jsp .* stor_mn.'));
        s2oc_odd(icp,idir) = -pi*trace((f.*delta_funct(o,ksi)) * imag(jorb .* stor_mn.'));
    end
end

end

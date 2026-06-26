function OptSpTork = OptSpintork_inj(H,r,o,f,p,wvf,nw,hv,direct,ksi)
% This function evaluates the optically driven spin torque
% The formula is quite similar as injection current, while not velocity
% difference, but spin torque difference. The spin torque is defined as
% T = [sigma, H]/i/hbar, where H is the Hamiltonian Hk
% One can also use Hsoc or Hxc to specify the soc or exchange torques

% Inputs:
% H - Hamiltonian, in size of nw x nw
% r - rmn, < m | r | n >
% o - energy difference, Em - En
% f - occupation difference, fm - fn
% p - projection matrix, in size of (nw x nw x nwl)
% wvf - wavefunction, in size of (nw x nw)
% nw - size of the Hamiltonian
% hv - optical frequency
% direct - the light polarization (i & j) and torque (a) directions (a, i, j)
% ksi - a broadening factor
nhv = length(hv);
ndir = size(direct,1);
nwl = size(p, 3);
% nlorb = zeros(ndir, nhv, nwl);
OptSpTork = zeros(ndir, nhv, nwl);
 
for idir = 1:ndir
a = direct(idir,1);
ele1 = direct(idir,2);
ele2 = direct(idir,3);
rb = r(:,:,ele1);
rc = r(:,:,ele2);
sptorka = calc_torque_sp(H,nw,wvf,a);

if ele1 == ele2  % LPL, zero for T and PT
    bc = real(rb.*rc.'+rc.*rb.');
else       % CPL, allowed for all general symmetry
    bc = -imag(rb.*rc.'-rc.*rb.');
end
% bc = real(rb.*rc.'+rc.*rb.');
for iwl = 1:nwl
sptorkp = (sptorka*p(:,:,iwl)+p(:,:,iwl)*sptorka)/2; % projection
diffsptork = diag(sptorkp)*ones(1,nw) - ones(nw,1)*diag(sptorkp).';
facsk1 = f.* real(diffsptork);
    
    for ihv = 1:nhv
        omega = hv(ihv);
        fac2 = bc .* delta_funct(o-omega,ksi);
        % nlorb(idir,ihv,iwl) = -pi*trace(facjo1*fac2)/2;
        OptSpTork(idir,ihv,iwl) = -pi*trace(facsk1*fac2)/2;
    end
end
end

end

function T = calc_torque_sp(H,nw,wavefunct,direct)
% We calculate the spin torque for each direction
% The torque is defined as T = 1/ihbar * [spin, H]
% Here H is the Hamiltonian. If it is just the exchange field, then you are
% ending up with the exchange torque; if it is the soc part, then just the
% soc torque. In general, we can give the total Hamiltonian, with the total
% torque effect
[sx,sy,sz] = build_Pauli(nw);
sx = sx/2;
sy = sy/2;
sz = sz/2;
Tx = comm(sx, H, -1)/1i;
Ty = comm(sy, H, -1)/1i;
Tz = comm(sz, H, -1)/1i;
if nargin == 3
    T = zeros(nw,nw,3);
    T(:,:,1) = wavefunct' * Tx * wavefunct;
    T(:,:,2) = wavefunct' * Ty * wavefunct;
    T(:,:,3) = wavefunct' * Tz * wavefunct;
elseif nargin == 4
    if direct == 1
        T = wavefunct' * Tx * wavefunct;
    elseif direct == 2
        T = wavefunct' * Ty * wavefunct;
    elseif direct == 3
        T = wavefunct' * Tz * wavefunct;
    else
        error('The direction of torque must be 1 - 3')
    end
else
    error('Number of input arguments wrongly given. Only 3 or 4 allowed')
end


end

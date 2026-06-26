% This code is read and analyze wannier90_hr.dat file and derive Floquet
% band structures
% Light is in vector potential form, A_i = A_0*cos(wt+phi_i),
% eta_i = e^(i*phi_i)
% left CPL: eta = [1, -i], A = [A0, A0]
% right CPL: eta = [1, i], A = [A0, A0]
% x-LPL: eta = [1, 0], A = [A0, 0]
% y-LPL: eta = [0, 1], A = [0, A0]
% Everything in eV or Angstrom
% Author: Jian Zhou
% Date: Mar 17, 2025
clear;

%% Floquet parameters
A0 = 0.2; % in eA/hbar [Ang^-1]
omega = 8; % in eV   (E[V/A] = 0.0005064 * A0 * omega)
eta = [1, -1i, 0];
% use eta = [1, -1i, 0] for left-CPL; eta = [1, 1i, 0] for right-CPL
% use eta = [1, 0, 0] for x-LPL; eta = [0, 1, 0] for y-LPL

nq = 1; % cutoff H_q (usually H^{-2} to H^{2})
A0 = A0 * eta;

%% Read Files
[R,HR,wtR,wc,wl,nw,nR] = getwanhr_layer('wannier90');
[alatt, blatt] = readlattice('POSCAR');
Vlatt = det(alatt);

%% Green parameters
kptpath = [0    0   0;
          -1/3 2/3  0;
          1/3  1/3  0;
           0   1/2  0;
           0   0    0];
KLABEL = ['G'; 'K';'K';'M';'G'];
ksep = 50; % how many points in each k-path section
EF = -2.0; % Fermi level

%% Main
[kptcell, kcorr, nkpts, klines] = k_path_highsym(kptpath, ksep, blatt);
Enk = zeros(nkpts,nw);
sz_k = zeros(nkpts, nw);
parpool(4)
parfor ikpt = 1:nkpts
    % disp(['ikpt: ',num2str(ikpt),' in total ',num2str(nkpts),' points'])
    kpt = kptcell{ikpt};
    [Hk,dHdk] = Magnus_Hk(R,HR,wtR,wc,nw,kpt,A0,omega,nq,alatt,blatt);
    % Hk = Magnus_Hk(R,HR,wtR,wc,nw,kpt,A0,omega,nq,alatt,blatt);
    [V, D] = eig(Hk);
    [di,ind] = sort(real(diag(D)));
    energy = real(di);
    wavefunct = V(:,ind);
    sz_loop = calc_sp(nw, wavefunct, 3);
    Enk(ikpt,:) = energy;
    sz_k(ikpt,:) = real(diag(sz_loop));
end
delete(gcp('nocreate'))

%% plotting
Eg = plotting_banddispersion(kcorr,Enk,EF,nw,nkpts,KLABEL,klines,ksep,sz_k);

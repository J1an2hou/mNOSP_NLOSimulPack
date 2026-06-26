% This code is read and analyze wannier90_hr.dat file
% Everything in eV or Angstrom
% We evaluate the polarization according to Berry phase method
clear;

%% Initial
nkx = 151;
nky = 151;
nkz = 1;
nocc = 20;
ksi = 0.02;

[R,HR,wtR,wc,nw,nR] = getwanhr('wannier90');
[alatt, blatt] = readlattice('POSCAR');
Vlatt = det(alatt);
Slatt = det(alatt(1:2,1:2));

%% Main
[kptcell, nkpts] = k_mesh_BZ(nkx, nky, nkz);
% [kptcell, nkpts] = generate_kmesh(nkx, nky, nkz);
bkx = zeros(1,nkpts);
bky = zeros(1,nkpts);
bkz = zeros(1,nkpts);
kweight = 1/nkpts;
f = zeros(nw, 1);
Enk = zeros(nkpts, nw);
wv_k = zeros(nkpts, nw, nw);

parpool(4)
parfor ikpt = 1:nkpts
    disp(['ikpt: ',num2str(ikpt),' in total ',num2str(nkpts),' points'])
    kpt = kptcell{ikpt};
    Hk = buildHk(R,HR,wtR,wc,nw,nR,kpt,alatt,blatt);
    [Vf, D_t] = eig(Hk);
    [d_sort,ind] = sort(real(diag(D_t)));
    energy = real(d_sort);
    wavefunct = Vf(:,ind);
    
    % k dependent quantities
    bkpt = kpt * blatt;
    Enk(ikpt, :) = energy;
    wv_k(ikpt, :, :) = wavefunct;
    bkx(ikpt) = bkpt(1);
    bky(ikpt) = bkpt(2);
    bkz(ikpt) = bkpt(3);
    
end
delete(gcp('nocreate'))
[Eg, EF] = band_gap_EF(Enk, nocc, 300);

Pdir1 = Polarization_Wilson(wv_k, nocc, nkx, nky, nkz, 1);
Pdir2 = Polarization_Wilson(wv_k, nocc, nkx, nky, nkz, 2);


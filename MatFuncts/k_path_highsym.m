function [kptcell, kcorr, kcnt, klines] = k_path_highsym(kptpath, ksep, blatt)
% Generate k coordinates along a high symmetric k-path
ks = size(kptpath,1);
klines = ks-1; % how many paths
nkpts = klines * ksep;
kpt_path = zeros(nkpts,3);
kptcell = cell(nkpts,1);
for iks = 1:klines
    kx = linspace(kptpath(iks,1),kptpath(iks+1,1),ksep);
    ky = linspace(kptpath(iks,2),kptpath(iks+1,2),ksep);
    kz = linspace(kptpath(iks,3),kptpath(iks+1,3),ksep);
    kpt_path((iks-1)*ksep+1:(iks-1)*ksep+ksep,1) = kx;
    kpt_path((iks-1)*ksep+1:(iks-1)*ksep+ksep,2) = ky;
    kpt_path((iks-1)*ksep+1:(iks-1)*ksep+ksep,3) = kz;
end
kcorr = zeros(1,nkpts);
for ikpt = 1:nkpts
    kpt = kpt_path(ikpt,:);
    bkpt = kpt*blatt;
    if ikpt == 1
        k_old = bkpt;
        kdisp = 0; % The horizontal axis coordinate of the k-path
        kmove = 0;
        kcnt = 1;
        kptcell{kcnt} = [kpt(1), kpt(2), kpt(3)];
    else
        kmove = bkpt - k_old;
        kmove = sqrt(kmove * kmove');
    end
    if ikpt ~= 1 && kmove < 1e-4
        continue
    elseif ikpt ~= 1
        kcnt = kcnt + 1;
        kdisp = kdisp + kmove;
        kcorr(kcnt) = kdisp;
        kptcell{kcnt} = [kpt(1), kpt(2), kpt(3)];
    end
    k_old = bkpt;
end
kptcell = kptcell(1:kcnt);
kcorr = kcorr(1:kcnt);

end
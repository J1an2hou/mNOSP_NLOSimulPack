function nlsp = calc_nle_inj(r,o,f,p,wvf,nw,hv,direct,ksi)
nhv = length(hv);
ndir = size(direct,1);
nwl = size(p, 3);
% nlorb = zeros(ndir, nhv, nwl);
nlsp = zeros(ndir, nhv, nwl);
 
for idir = 1:ndir
a = direct(idir,1);
ele1 = direct(idir,2);
ele2 = direct(idir,3);
rb = r(:,:,ele1);
rc = r(:,:,ele2);
% oa = calc_orb(nw,wvf,orbital_seq,a);
sa = calc_sp(nw,wvf,a);
if ele1 == ele2  % LPL, zero for T and PT
    bc = real(rb.*rc.'+rc.*rb.');
else       % CPL, allowed for all general symmetry
    bc = -imag(rb.*rc.'-rc.*rb.');
end
for iwl = 1:nwl
% oap = (oa*p(:,:,iwl)+p(:,:,iwl)*oa)/2; % projection
sap = (sa*p(:,:,iwl)+p(:,:,iwl)*sa)/2; % projection
% diffod = diag(oap)*ones(1,nw) - ones(nw,1)*diag(oap).';
diffsd = diag(sap)*ones(1,nw) - ones(nw,1)*diag(sap).';
% facjo1 = f.* real(diffod);
facjs1 = f.* real(diffsd);
    
    for ihv = 1:nhv
        omega = hv(ihv);
        fac2 = bc .* delta_funct(o-omega,ksi);
        % nlorb(idir,ihv,iwl) = -pi*trace(facjo1*fac2)/2;
        nlsp(idir,ihv,iwl) = -pi*trace(facjs1*fac2)/2;
    end
end
end

end


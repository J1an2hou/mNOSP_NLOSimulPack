function sp_tork_shift = OptSpintork_shift(o,r,f,p,wvf,nw,hv,ksi,direct)
nhv = length(hv);
ndir = size(direct,1);
nwl = size(p,3);
sp_tork_shift = zeros(ndir, nhv, nwl);

for idir = 1:ndir
a = direct(idir,1);
ele1 = direct(idir,2);
ele2 = direct(idir,3);
sa = calc_sp(nw,wvf,a);
% sptorka = calc_torque_sp(H,nw,wvf,a);
for iwl = 1:nwl
rt_abc = calc_rs(r,p(:,:,iwl),ele1,ele2,sa);
rt_acb = calc_rs(r,p(:,:,iwl),ele2,ele1,sa);
if ele1 == ele2   % LPL
    facs1 = -imag(rt_abc-rt_acb.');
else       % CPL 
    facs1 = real(rt_abc+rt_acb.');
end
for ihv = 1:nhv
    omega = hv(ihv);
    fac2 = delta_funct(o-omega,ksi);
    % if ele1 == ele2 
    %     fac2 = delta_funct(o-omega,ksi)+delta_funct(o+omega,ksi);
    % else    
    %     fac2 = delta_funct(o-omega,ksi)-delta_funct(o+omega,ksi);
    % end
    sp_tork_shift(idir,ihv,iwl) = pi*trace((f.*fac2)*facs1)/4;
end
end
end
end

function dsdk = calc_dsdk(r,p,c,sa)
% sigma^a_;c = i*(r_c*sigma_a - sigma_a*r_c)
rc = r(:,:,c);
sa = (sa*p+p*sa)/2;

dsdk = 1i*(rc * sa - sa * rc);
end

function rs_bcsa = calc_rs(r,p,b,c,sa)
% index (sa,b,c): r_mn^b * rcsa_nm
rb = r(:,:,b);
rc_sa = calc_dsdk(r,p,c,sa);
rs_bcsa = rb.' .* rc_sa;
end
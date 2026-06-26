function nlsp = calc_nle_shift(o,v,r,f,p,wft,nw,hv,ksi,direct)
nhv = length(hv);
ndir = size(direct,1);
nwl = size(p,3);
nlsp = zeros(ndir, nhv, nwl);
% nlorb = zeros(ndir, nhv, nwl);

for idir = 1:ndir
a = direct(idir,1);
ele1 = direct(idir,2);
ele2 = direct(idir,3);
sa = calc_sp_updn(nw,wft,a);
% oa = calc_orb(nw,wft,atomorb,a);
for iwl = 1:nwl
rs_abc = calc_rs(r,v,o,p(:,:,iwl),ksi,ele1,ele2,sa);
rs_acb = calc_rs(r,v,o,p(:,:,iwl),ksi,ele2,ele1,sa);
% ro_abc = calc_rs(r,v,o,p(:,:,iwl),ksi,ele1,ele2,oa);
% ro_acb = calc_rs(r,v,o,p(:,:,iwl),ksi,ele2,ele1,oa);
if ele1 == ele2   % Normal Shift NLE (LPL, allowed for P, T, PT)
    facs1 = imag(rs_abc+rs_acb);
    % faco1 = imag(ro_abc+ro_acb);
else       % Magnetic Shift NLE (CPL, zero for T and for PT)
    facs1 = real(rs_abc-rs_acb);
    % faco1 = real(ro_abc-ro_acb);
end
for ihv = 1:nhv
    omega = hv(ihv);
    if ele1 == ele2   % Normal Shift NLE
        fac2 = delta_funct(o-omega,ksi)+delta_funct(o+omega,ksi);
    else       % Magnetic Shift NLE
        fac2 = delta_funct(o-omega,ksi)-delta_funct(o+omega,ksi);
    end
    nlsp(idir,ihv,iwl) = pi*trace((f.*fac2)*facs1)/4;
    % nlorb(idir,ihv,iwl) = pi*trace((f.*fac2)*faco1)/4;
end
end
end
end

function rdk = calc_dsdk(v,o,p,eta,c,sa)
% rcsd_nm^{c;(sd)}
vc = v(:,:,c);
sa = (sa*p+p*sa)/2;
% Dc = diag(real(vc))*ones(1,nw)-ones(nw,1)*diag(real(vc))';
% rdk = 1i*((sd.*Dc).*o./(o.^2+eta^2) ...
%     +(vc*(sd.*o./(o.^2+eta^2)) - (sd.*o./(o.^2+eta^2))*vc)) ...
%     .* o./(o.^2+eta^2);
rdk = 1i*(vc*(sa.*o./(o.^2+eta^2)) - (sa.*o./(o.^2+eta^2))*vc) ...
    .* o./(o.^2+eta^2);
end

function rs_bcsa = calc_rs(r,v,o,p,eta,b,c,sa)
% index (sa,b,c): r_mn^b * rcsd_nm^{c,sa}
rb = r(:,:,b);
rc_sa = calc_dsdk(v,o,p,eta,c,sa);
rs_bcsa = rb.' .* rc_sa;
end

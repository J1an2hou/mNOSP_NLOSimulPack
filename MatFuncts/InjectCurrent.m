function [bopv_on_k,bspv_on_k] ...
    = InjectCurrent(v,r,o,f,p,wvf,orbital_seq,nw,hv,direct,ksi)
nhv = length(hv);
ndir = size(direct,1);
nwl = size(p, 3);
bopv_on_k = zeros(ndir, nhv, nwl);
bspv_on_k = zeros(ndir, nhv, nwl);
%jorb = zeros(ndir,nw);
%jspin = zeros(ndir,nw);
 
for idir = 1:ndir
    a = direct(idir,1);
    b = direct(idir,2);
    c = direct(idir,3);
    d = direct(idir,4);
va = v(:,:,a);
ra = r(:,:,a);
od = calc_orb(nw,wvf,orbital_seq,d);
sd = calc_sp(nw,wvf,d);
jod = (od*va+va*od)/2;
jod = jod + 1i*(ra*(o.*od)+(o.*od)*ra)/2; % This is the added torque term
jsd = (sd*va+va*sd)/2;
jsd = jsd + 1i*(ra*(o.*sd)+(o.*sd)*ra)/2; % This is the added torque term
if b == c  % Magnetic Injection Current
          %(LPL, charge current: zero for T; spin current: zero for PT)
    bc = real(r(:,:,b).*r(:,:,c).'+r(:,:,c).*r(:,:,b).');
else       % Normal Injection Current
          %(CPL, charge current: zero for PT; spin current: zero for T)
    bc = imag(r(:,:,b).*r(:,:,c).'-r(:,:,c).*r(:,:,b).');
end
for iwl = 1:nwl
jodp = jod*p(:,:,iwl); % projection
jsdp = jsd*p(:,:,iwl); % projection
diffjod = diag(jodp)*ones(1,nw) - ones(nw,1)*diag(jodp).';
diffjsd = diag(jsdp)*ones(1,nw) - ones(nw,1)*diag(jsdp).';
%jorb(idir,:) = diag(real(jod));
%jspin(idir,:) = diag(real(jsd));
 
facjo1 = f.* real(diffjod);
facjs1 = f.* real(diffjsd);
    
    for ihv = 1:nhv
        omega = hv(ihv);
        fac2 = bc .* df(o-omega,ksi);
        bopv_on_k(idir,ihv,iwl) = -trace(fac2*facjo1)/2;
        bspv_on_k(idir,ihv,iwl) = -trace(fac2*facjs1)/2;
    end
end
end

end
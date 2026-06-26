function [SC_sk,SC_ok,Rnmab_s,Rnmab_o] ...
    = ShiftCurrent_SurfaceZ(o,v,w,r,f,p,wavefunct,orbital_seq,nw,nL,hv,ksi,direct)
nhv = length(hv);
ndir = size(direct,1);
nwl = size(p, 3);
SC_sk = zeros(ndir, nhv,nwl);
SC_ok = zeros(ndir, nhv,nwl);
Rnmab_s = zeros(nw*nL,nw*nL,ndir,nwl);  % shift vector
Rnmab_o = zeros(nw*nL,nw*nL,ndir,nwl);
for idir = 1:ndir
a = direct(idir,1);
b = direct(idir,2);
c = direct(idir,3);
d = direct(idir,4);
sd = calc_sp(nw,nL,wavefunct,d);
od = calc_orb(nw,nL,wavefunct,orbital_seq,d);
for iwl = 1:nwl
Iabcsd = calc_Imn(r,v,o,w,p(:,:,iwl),nw*nL,ksi,a,b,c,sd);
Iacbsd = calc_Imn(r,v,o,w,p(:,:,iwl),nw*nL,ksi,a,c,b,sd);
Iabcod = calc_Imn(r,v,o,w,p(:,:,iwl),nw*nL,ksi,a,b,c,od);
Iacbod = calc_Imn(r,v,o,w,p(:,:,iwl),nw*nL,ksi,a,c,b,od);
if b == c   % Normal Shift Current (LPL, zero for PT)
    facs1 = imag(Iabcsd+Iacbsd);
    faco1 = imag(Iabcod+Iacbod);
    r2 = abs(r(:,:,b)).^2;
    Rnmab_s(:,:,idir,iwl) = imag(Iabcsd) .* r2 ./ (r2.^2 + ksi^2);
    Rnmab_o(:,:,idir,iwl) = imag(Iabcod) .* r2 ./ (r2.^2 + ksi^2);
else       % Magnetic Shift Current (CPL, zero for T)
    facs1 = real(Iabcsd-Iacbsd);
    faco1 = real(Iabcod-Iacbod);
end

for ihv = 1:nhv
    omega = hv(ihv);
    if b == c   % Normal Shift Current (LPL, zero for PT)
        fac2 = df(o-omega,ksi)+df(o+omega,ksi);
    else       % Magnetic Shift Current (CPL, zero for T)
        fac2 = df(o-omega,ksi)-df(o+omega,ksi);
    end
    SC_sk(idir,ihv,iwl) = trace((f.*fac2)*facs1)/4;
    SC_ok(idir,ihv,iwl) = trace((f.*fac2)*faco1)/4;
end
end

end
end
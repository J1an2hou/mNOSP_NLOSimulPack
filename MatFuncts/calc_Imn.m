function Iabcsd = calc_Imn(r,v,o,w,p,nw,eta,a,b,c,sd)
% for shift current calculation
% index (a,b,c,sd): r_mn^b * r_nm^{c;(a,sd)}
rb = r(:,:,b);
rc_dkad = calc_rdk(v,o,w,p,nw,eta,c,a,sd);
Iabcsd = rb.' .* rc_dkad;
end
function rdk = calc_rdk(v,o,w,p,nw,eta,a,b,sd)
% Covarent gradient of position over k
% Sum rule method is applied
% r_nm^{a;(b,spind)} or d r_nm^{a} / d k_{b,spind}

if nargin == 8 % without spin index clearly indicated
    sd = eye(nw, nw);
end

va = v(:,:,a);
vb = v(:,:,b);
vb = (vb*sd+sd*vb)/2*p;
wab = w(:,:,a,b);
wab = (wab*sd+sd*wab)/2*p;
Da = diag(real(va))*ones(1,nw)-ones(nw,1)*diag(real(va))';
% Db = diag(real(vb))*ones(1,nw)-ones(nw,1)*diag(real(vb))';
% rdk = 1i*((va.*Db+vb.*Da).*o./(o.^2+eta^2) ...
%     +(va*(vb.*o./(o.^2+eta^2)) - (vb.*o./(o.^2+eta^2))*va) ...
%     -wab) .* o./(o.^2+eta^2);
rdk = 1i*((vb.*Da).*o./(o.^2+eta^2) ...
    +(va*(vb.*o./(o.^2+eta^2)) - (vb.*o./(o.^2+eta^2))*va) ...
    -wab) .* o./(o.^2+eta^2);
end
function [r, o] = calc_rmn_omn(v,energy,nw,eta)
r = zeros(nw,nw,3);
o = energy*ones(1,nw) - ones(nw,1)*energy';
v1 = v(:,:,1);
v2 = v(:,:,2);
v3 = v(:,:,3);
r1 = -1i*v1 .* o ./ (o.^2 + eta^2);
r2 = -1i*v2 .* o ./ (o.^2 + eta^2);
r3 = -1i*v3 .* o ./ (o.^2 + eta^2);
r(:,:,1) = r1;
r(:,:,2) = r2;
r(:,:,3) = r3;
end
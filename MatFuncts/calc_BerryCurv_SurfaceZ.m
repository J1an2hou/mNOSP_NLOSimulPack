function [bc,bc_on_k] = calc_BerryCurv_SurfaceZ(v, o, p, f, wavefunct, nw, nL, ksi)
% We calculate the Spin Berry curvature Im(sp*vx, vy), here sp is
% 4-dimensional vector. Thus, the bc is also 4-dimensional
nwl = size(p, 3);
bc = zeros(4, nwl);
bc_on_k = zeros(nw*nL,4,nwl);
vx = v(:,:,1);
vy = v(:,:,2);
ry = vy .* o ./(o.^2+ksi^2);
for i = 0:3
spin = calc_sp(nw,nL,wavefunct,i);
vxsp = comm(vx,spin,+1)/2;
rxsp = vxsp .* o ./(o.^2+ksi^2);
for iwl = 1:nwl
rxsp_projl = comm(rxsp,p(:,:,iwl),+1)/2;
bc(i+1, iwl) = trace(imag((f.*ry) * rxsp_projl));
bc_on_k(:,i+1,iwl) = diag(2*imag(rxsp_projl*ry));
end
end
 
end
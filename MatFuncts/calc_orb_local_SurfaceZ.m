function o = calc_orb_local_SurfaceZ(nw,nL,wavefunct,orbital_seq,direct)
L = eye(nL);
[Jx, Jy, Jz] = buildJ(orbital_seq);
%Jx = projAMonbands(Jx,wavefunct,nw,nL);
%Jy = projAMonbands(Jy,wavefunct,nw,nL);
%Jz = projAMonbands(Jz,wavefunct,nw,nL);
if direct == 0
    o = wavefunct'*kron(L, eye(nw))*wavefunct;
elseif direct == 1
    o = wavefunct'*kron(L, Jx)*wavefunct;
elseif direct == 2
    o = wavefunct'*kron(L, Jy)*wavefunct;
elseif direct == 3
    o = wavefunct'*kron(L, Jz)*wavefunct;
end
 
end
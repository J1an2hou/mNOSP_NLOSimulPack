function o = calc_orb_local(nw,wavefunct,orbital_seq,direct)
A = eye(nw,nw);
[Jx, Jy, Jz] = buildJ(orbital_seq);
%Jx = projAMonbands(Jx,wavefunct,nw,nL);
%Jy = projAMonbands(Jy,wavefunct,nw,nL);
%Jz = projAMonbands(Jz,wavefunct,nw,nL);
if nargin < 4
    o = zeros(nw,nw,3);
    o(:,:,1) = wavefunct' * Jx * wavefunct;
    o(:,:,2) = wavefunct' * Jy * wavefunct;
    o(:,:,3) = wavefunct' * Jz * wavefunct;
else

    if direct == 0
    o = wavefunct'* A *wavefunct;
    elseif direct == 1
    o = wavefunct'* Jx *wavefunct;
    elseif direct == 2
    o = wavefunct'* Jy *wavefunct;
    elseif direct == 3
    o = wavefunct'* Jz *wavefunct;
    end
end
 
end
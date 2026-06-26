function s = calc_sp(nw,wavefunct,direct)
[sx,sy,sz] = build_Pauli(nw);
sx = sx/2;
sy = sy/2;
sz = sz/2;
if nargin == 2
    s = zeros(nw,nw,3);
    s(:,:,1) = wavefunct' * sx * wavefunct;
    s(:,:,2) = wavefunct' * sy * wavefunct;
    s(:,:,3) = wavefunct' * sz * wavefunct;
elseif nargin == 3
    if direct == 0
        s = wavefunct' * wavefunct;
    elseif direct == 1
        s = wavefunct' * sx * wavefunct;
    elseif direct == 2
        s = wavefunct' * sy * wavefunct;
    elseif direct == 3
        s = wavefunct' * sz * wavefunct;
    else
        error('Spin direction must be 0 - 3')
    end
else
    error('Number of input arguments wrongly given. Only 2 or 3 allowed')
end


end

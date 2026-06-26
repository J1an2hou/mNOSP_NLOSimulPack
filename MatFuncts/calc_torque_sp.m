function T = calc_torque_sp(H,nw,wavefunct,direct)
% We calculate the spin torque for each direction
% The torque is defined as T = 1/ihbar * [spin, H]
% Here H is the Hamiltonian. If it is just the exchange field, then you are
% ending up with the exchange torque; if it is the soc part, then just the
% soc torque. In general, we can give the total Hamiltonian, with the total
% torque effect
[sx,sy,sz] = build_Pauli(nw);
sx = sx/2;
sy = sy/2;
sz = sz/2;
Tx = comm(sx, H, -1)/1i;
Ty = comm(sy, H, -1)/1i;
Tz = comm(sz, H, -1)/1i;
if nargin == 3
    T = zeros(nw,nw,3);
    T(:,:,1) = wavefunct' * Tx * wavefunct;
    T(:,:,2) = wavefunct' * Ty * wavefunct;
    T(:,:,3) = wavefunct' * Tz * wavefunct;
elseif nargin == 4
    if direct == 1
        T = wavefunct' * Tx * wavefunct;
    elseif direct == 2
        T = wavefunct' * Ty * wavefunct;
    elseif direct == 3
        T = wavefunct' * Tz * wavefunct;
    else
        error('The direction of torque must be 1 - 3')
    end
else
    error('Number of input arguments wrongly given. Only 3 or 4 allowed')
end


end

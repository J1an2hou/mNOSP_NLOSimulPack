
function s = Pauli_mat(nw, dir)
A = eye(nw/2);
sigmax = [0 1; 1 0];
sigmay = [0 -1i; 1i 0];
sigmaz = [1 0; 0 -1];
sigma0 = [1 0; 0 1];
% The basis function is |1,up>, |2,up>, ..., |1,dn>, |2,dn>, ...
s0 = kron(sigma0, A);
sx = kron(sigmax, A);
sy = kron(sigmay, A);
sz = kron(sigmaz, A);
if dir == 0
    s = s0;
elseif dir == 1
    s = sx;
elseif dir == 2
    s = sy;
elseif dir == 3
    s = sz;
else
    error('Error building Pauli matrix')
end

end


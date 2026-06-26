function s = calc_sp_surfaceZ_remove(nw,nL,rmindx,wavefunct,direct)
s = zeros(nw,nw);
A = eye(nw/2);
L = eye(nL);
sigmax = [0 1; 1 0];
sigmay = [0 -1i; 1i 0];
sigmaz = [1 0; 0 -1];
sigma0 = [1 0; 0 1];
if direct == 0
    s = kron(L, kron(sigma0, A));
elseif direct == 1
    s = kron(L, kron(sigmax, A))/2;
elseif direct == 2
    s = kron(L, kron(sigmay, A))/2;
elseif direct == 3
    s = kron(L, kron(sigmaz, A))/2;
end
s = MatTrim_RowsCols(s,rmindx);
s = wavefunct'*s*wavefunct;
end
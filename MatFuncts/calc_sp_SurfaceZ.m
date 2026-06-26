function s = calc_sp_SurfaceZ(nw,nL,wavefunct,direct)
% s = zeros(nw,nw);
A = eye(nw/2);
L = eye(nL);
sigmax = [0 1; 1 0]/2;
sigmay = [0 -1i; 1i 0]/2;
sigmaz = [1 0; 0 -1]/2;
sigma0 = [1 0; 0 1];
if nargin == 3 % no specified direction indicated
    s = zeros(nw*nL,nw*nL,3);
    s(:,:,1) = wavefunct' * kron(L, kron(sigmax, A)) * wavefunct;
    s(:,:,2) = wavefunct' * kron(L, kron(sigmay, A)) * wavefunct;
    s(:,:,3) = wavefunct' * kron(L, kron(sigmaz, A)) * wavefunct;
elseif nargin == 4
    if direct == 0
        s = kron(L, kron(sigma0, A));
    elseif direct == 1
        s = kron(L, kron(sigmax, A));
    elseif direct == 2
        s = kron(L, kron(sigmay, A));
    elseif direct == 3
        s = kron(L, kron(sigmaz, A));
    else
        error('Spin direction must be 0 - 3')
    end
else
    error('Number of input arguments wrongly given. Only 3 or 4 allowed')
end

end
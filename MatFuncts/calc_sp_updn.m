function s = calc_sp_updn(nw,wavefunct,direct)
s = zeros(nw,nw);
A = eye(nw/2);
sigmax = [0 1; 1 0];
sigmay = [0 -1i; 1i 0];
sigmaz = [1 0; 0 -1];
sigma0 = [1 0; 0 1];
% The basis function is |1,up>, |2,up>, ..., |1,dn>, |2,dn>, ...
% if direct == 0
%     s = kron(sigma0, A);
% elseif direct == 1
%     s = kron(sigmax, A)/2;
% elseif direct == 2
%     s = kron(sigmay, A)/2;
% elseif direct == 3
%     s = kron(sigmaz, A)/2;
% end

% The basis function is |1,up>, |1,dn>, |2,up>, |2,dn>, ...
if direct == 0
    s = kron(A, sigma0);
elseif direct == 1
    s = kron(A, sigmax)/2;
elseif direct == 2
    s = kron(A, sigmay)/2;
elseif direct == 3
    s = kron(A, sigmaz)/2;
end


s = wavefunct'*s*wavefunct;

end

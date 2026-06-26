function s = calc_neel(wavefunct,site_updn,direct)

A_size = size(site_updn,1);
Nw = size(wavefunct,1);
if Nw == A_size
    A = site_updn(1:2:Nw,1:2:Nw);
elseif Nw == A_size*2
    A = site_updn;
else
    error('Wrong size of site_updn matrix')
end
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


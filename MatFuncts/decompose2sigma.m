function [H0,H1,H2,H3] = decompose2sigma(Hk,nw, nL)
% Hk is in size of (nw x nw) with basis set |orb1, up>, |orb2, up> ... to
% |orb_n, up>, then |orb_1, dn>, |orb_2, dn> ... |orb_n, dn>
% This function is to decompose a Hermitian Hk according to Pauli matrices
% The parameter nL is optional. If it is larger than 1, then we duplicate
% the basis set according to each layer (nL). In that case, the orbital sequence is
% (|1,up>,|2,up>,|1,dn>,|2,dn>)@layer-1, then
% (|1,up>,|2,up>,|1,dn>,|2,dn>)@layer-2, ..., until layer-nL

if nargin == 2
    nL = 1;
end
L = eye(nL);
rnd = -6;
Hk = roundn(Hk, rnd);
if ~ishermitian(Hk)
    warning('The input Hamiltonian is not Hermitian, danger to do decomposition')
end
nhalf = nw/2;
sx = [0, 1; 1, 0];
sy = [0, -1i; 1i, 0];
sz = [1, 0; 0, -1];
% s0 = [1, 0; 0, 1];
sigmax = kron(L, kron(sx, eye(nhalf, nhalf)));
sigmay = kron(L, kron(sy, eye(nhalf, nhalf)));
sigmaz = kron(L, kron(sz, eye(nhalf, nhalf)));
% sigma0 = kron(s0, eye(nhalf, nhalf));
% Now let's do some mapping according to sigmax * Hk * sigmax = H0 + H1 ...
xHx = sigmax * Hk * sigmax;
g01 = (Hk + xHx)/2;
g23 = (Hk - xHx)/2;
yg01y = sigmay * g01 * sigmay;
zg23z = sigmaz * g23 * sigmaz;
H0 = (g01 + yg01y)/2;
H1 = (g01 - yg01y)/2;
H3 = (g23 + zg23z)/2;
H2 = (g23 - zg23z)/2;

H0 = roundn(H0, rnd);
H1 = roundn(H1, rnd);
H2 = roundn(H2, rnd);
H3 = roundn(H3, rnd);
% Hdiff = Hk - H0 - H1 - H2 - H3;
% if norm(Hdiff) > 1e-6
%     warning(['The decomposition of Hamiltonian may be incomplete: ',num2str(norm(Hdiff)),' is the residual'])
% end

end
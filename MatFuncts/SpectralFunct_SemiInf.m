function [SurfA0,SurfAs] = SpectralFunct_SemiInf(H00,H01,epsilon,layer,nw,nL,ksi)
neps = length(epsilon);
nlayer = length(layer);
SurfA0 = zeros(neps,nlayer); % Top Strip spectral function
SurfAs = zeros(neps,3,nlayer); % Strip spectral function for spin-1,2,3
% N = nw*nL;
% Grnn_top = zeros(N,N,nlayer,neps);
sp1 = sp4Spectral(nw,nL,1);
sp2 = sp4Spectral(nw,nL,2);
sp3 = sp4Spectral(nw,nL,3);
for ieps = 1:neps
    eps = epsilon(ieps);
    Gtop = SurfGreen(H00,H01,eps+1i*ksi); % retarded Green function
    g0 = InnerGreen(Gtop,H01,layer);      % g0: N*N*nlayer
    % Grnn_top(:,:,:,ieps) = g0;
    for ilayer = 1:nlayer        
        SurfA0(ieps,ilayer) = -imag(trace(g0(:,:,ilayer)))/pi;
        % SurfAs(ieps,1,ilayer) = -imag(trace(sp1*g0(:,:,ilayer))) ...
        %     *SurfA0(ieps,ilayer)/(SurfA0(ieps,ilayer)^2+ksi^2)/pi;
        % SurfAs(ieps,2,ilayer) = -imag(trace(sp2*g0(:,:,ilayer))) ...
        %     *SurfA0(ieps,ilayer)/(SurfA0(ieps,ilayer)^2+ksi^2)/pi;
        % SurfAs(ieps,3,ilayer) = -imag(trace(sp3*g0(:,:,ilayer))) ...
        %     *SurfA0(ieps,ilayer)/(SurfA0(ieps,ilayer)^2+ksi^2)/pi;
        SurfAs(ieps,1,ilayer) = -imag(trace(sp1*g0(:,:,ilayer)));
        SurfAs(ieps,2,ilayer) = -imag(trace(sp2*g0(:,:,ilayer)));
        SurfAs(ieps,3,ilayer) = -imag(trace(sp3*g0(:,:,ilayer)));
    end
end

end

function s = sp4Spectral(nw,nL,direct)
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
end



function Gtop = SurfGreen(H00,H01,E)
maxcyc = 100;
tol = 1e-10;
eps = H00; epstop = H00; epsbot = H00; alpha = H01; beta = H01';
ndim = size(H00,1);
omega = E*eye(ndim,ndim);
for icyc = 1:maxcyc
    alpha_new = alpha/(omega-eps)*alpha;
    beta_new = beta/(omega-eps)*beta;
    eps_new = eps+alpha/(omega-eps)*beta+beta/(omega-eps)*alpha;
    epstop_new = epstop+alpha/(omega-eps)*beta;
    epsbot_new = epsbot+beta/(omega-eps)*alpha;
    if norm(alpha_new,'fro') < tol
        break
    end
    alpha = alpha_new;
    beta = beta_new;
    eps = eps_new;
    epstop = epstop_new;
    epsbot = epsbot_new;
    if icyc == maxcyc
        disp('Max cycle of SurfGreen function reached. Warning!')
        disp(['The energy diff in bulk is still: ',num2str(norm(alpha_new,'fro'))])
    end
end
%Gbulk = inv(omega-eps);
Gtop = inv(omega-epstop);
%Gbot = inv(omega-epsbot);

end

function Gnntop = InnerGreen(Gt,H01,layer)
nlayer = length(layer);
layer = sort(layer);
ndim = size(Gt,1);
Gnntop = zeros(ndim,ndim,nlayer);
maxlayer = max(layer);
ifind = 1;
H10 = H01';
Gnn = Gt;
for il = 0:maxlayer
    if ismember(il,layer)
        Gnntop(:,:,ifind) = Gnn;
        ifind = ifind+1;
    end
    Gnn = Gt+Gt*H10*Gnn*H01*Gt;
end

end
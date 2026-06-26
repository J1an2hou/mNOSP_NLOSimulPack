function SurfA0 = SpectralFunct_SemiInf_noSpin(H00,H01,epsilon,layer,ksi)
neps = length(epsilon);
nlayer = length(layer);
SurfA0 = zeros(neps,nlayer); % Top Strip spectral function
% Grnn_top = zeros(N,N,nlayer,neps);

for ieps = 1:neps
    eps = epsilon(ieps);
    Gtop = SurfGreen(H00,H01,eps+1i*ksi); % retarded Green function
    g0 = InnerGreen(Gtop,H01,layer);      % g0: N*N*nlayer
    % Grnn_top(:,:,:,ieps) = g0;
    for ilayer = 1:nlayer        
        SurfA0(ieps,ilayer) = -imag(trace(g0(:,:,ilayer)))/pi;
    end
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
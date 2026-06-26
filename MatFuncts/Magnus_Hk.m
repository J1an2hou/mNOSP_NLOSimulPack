function [Hk, dHdk] = Magnus_Hk(R,HR,wtR,wc,nw,kpt,A0,omega,nq,alatt,blatt)
% We build the Floquet Hk for a slab terminated along z direction
Hkq = zeros(nw,nw,nq);
Hk_q = zeros(nw,nw,nq);
dHdkq = zeros(nw,nw,3,nq);
dHdk_q = zeros(nw,nw,3,nq);
nR = size(R,1);
FloPhase0 = FloPhase_ave(A0,R,wc,alatt,nw,nR,0);
[Hkq0,dHdkq0] = buildHk(R,HR.*FloPhase0,wtR,wc,nw,nR,kpt,alatt,blatt);
for iq = 1:nq
    FloPhaseq = FloPhase_ave(A0,R,wc,alatt,nw,nR,iq);
    FloPhase_q = FloPhase_ave(A0,R,wc,alatt,nw,nR,-iq);
    [Hkq(:,:,iq),dHdkq(:,:,:,iq)] ...
        = buildHk(R,HR.*FloPhaseq,wtR,wc,nw,nR,kpt,alatt,blatt);
    [Hk_q(:,:,iq),dHdk_q(:,:,:,iq)] ...
        = buildHk(R,HR.*FloPhase_q,wtR,wc,nw,nR,kpt,alatt,blatt);
end
Hk = Hkq0;
dHdk = dHdkq0;
for iq = 1:nq
    Hk = Hk + comm(Hk_q(:,:,iq),Hkq(:,:,iq),-1)/iq/omega;
    if nargout > 1
    for i = 1:3
        dHdk(:,:,i) = dHdk(:,:,i) ...
            + comm(dHdk_q(:,:,i,iq),Hkq(:,:,iq),-1)/omega/iq ...
            + comm(Hk_q(:,:,iq),dHdkq(:,:,i,iq),-1)/omega/iq;
    end
    end
end
rdn = -8;
Hk = roundn(Hk,rdn);
if ~ishermitian(Hk)
    warning('Hamiltonian is not Hermitian!')
end
if nargout > 1
for i = 1:3
    dHdk(:,:,i) = roundn(dHdk(:,:,i), rdn);
    if ~ishermitian(dHdk(:,:,i))
        warning('Hamiltonian derivative is not Hermitian!')
    end
end
end
end
function Hk = Magnus_Hk_StripeX(R,HR,wtR,wc,nw,kpt,nL,A0,omega,nq,alatt,blatt)
Hkq = zeros(nw*nL,nw*nL,nq);
Hk_q = zeros(nw*nL,nw*nL,nq);
nR = size(R,1);
FloPhase0 = FloPhase_ave(A0,R,wc,alatt,nw,nR,0);
[Hkq0,~,~] = buildHk_StripeX(HR.*FloPhase0,wtR,R,wc,nw,kpt,alatt,blatt,nL);
for iq = 1:nq
    FloPhaseq = FloPhase_ave(A0,R,wc,alatt,nw,nR,iq);
    FloPhase_q = FloPhase_ave(A0,R,wc,alatt,nw,nR,-iq);
    [Hkq(:,:,iq),~,~] ...
        = buildHk_StripeX(HR.*FloPhaseq,wtR,R,wc,nw,kpt,alatt,blatt,nL);
    [Hk_q(:,:,iq),~,~] ...
        = buildHk_StripeX(HR.*FloPhase_q,wtR,R,wc,nw,kpt,alatt,blatt,nL);
end
Hk = Hkq0;
for iq = 1:nq
    Hk = Hk + comm(Hk_q(:,:,iq),Hkq(:,:,iq),-1)/iq/omega;
end
rdn = -8;
Hk = roundn(Hk,rdn);
if ~ishermitian(Hk)
    warning('Hamiltonian is not Hermitian!')
end
end
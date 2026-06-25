function [H00,H01] = BlockHam_SurfZ_FloquetA(HR,wtR,R,omega,A0,nq,wc,nw,kpt,alatt,blatt,nL)
% We build a supercell of nL unit cells, terminated along b-lattice
% This can be used to calculate surface state of a bullk system
% The orbital sequences are 1, ..., nw (in cell-1), 1, ..., nw (in cell-2)
% ..., until cell-nL.
% This nL unit cells form a principal layer for Green's function
nR = size(R,1);
R3max = max(R(:,3));
R3min = min(R(:,3));
nR3 = R3max-R3min+1;
HkqPn = zeros(nw,nw,nR3,nq);
Hk_qPn = zeros(nw,nw,nR3,nq);

FloPhase0 = FloPhase_ave(A0,R,wc,alatt,nw,nR,0);
[Hk_q0Pn,~,~]=buildHkAlongb1b2(R,HR.*FloPhase0,wtR,wc,nw,kpt,alatt,blatt);
% Hk_q0Pn is in size(nw,nw,nR2), same as dHdk_q0Pn (it is only dkx)
for iq = 1:nq
    FloPhaseq = FloPhase_ave(A0,R,wc,alatt,nw,nR,iq);
    FloPhase_q = FloPhase_ave(A0,R,wc,alatt,nw,nR,-iq);
    [HkqPn(:,:,:,iq),~,~] = buildHkAlongb1b2(R,HR.*FloPhaseq,wtR,wc,nw,kpt,alatt,blatt);
    [Hk_qPn(:,:,:,iq),~,~] = buildHkAlongb1b2(R,HR.*FloPhase_q,wtR,wc,nw,kpt,alatt,blatt);
end
HkPn = Hk_q0Pn;
for iq = 1:nq
    for iR3 = 1:nR3
    HkPn(:,:,iR3) = HkPn(:,:,iR3) ...
        + comm(Hk_qPn(:,:,iR3,iq),HkqPn(:,:,iR3,iq),-1)/omega/iq;
    end
end

H00 = zeros(nL*nw,nL*nw);
% V = zeros(nL*nw,nL*nw,3);
H01 = zeros(nL*nw,nL*nw);
R3max = max(R(:,3));
R3min = min(R(:,3));
for iL = 1:nL
for jL = 1:nL
    iR3 = jL-iL+1-R3min;
    iR3s = jL-iL+1-R3min-nL;
    if iR3 >=1 && iR3 <= R3max-R3min+1
        H00((iL-1)*nw+1:iL*nw,(jL-1)*nw+1:jL*nw) = HkPn(:,:,iR3);
        % V((iL-1)*nw+1:iL*nw,(jL-1)*nw+1:jL*nw,1) = dHdkPn(:,:,iR3,1);
    end
    if iR3s >= 1 && iR3s <= R3max-R3min+1
        H01((iL-1)*nw+1:iL*nw,(jL-1)*nw+1:jL*nw) = HkPn(:,:,iR3s);
    end
end
end

end
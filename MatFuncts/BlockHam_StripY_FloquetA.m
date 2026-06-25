
function [H00,H01] = BlockHam_StripY_FloquetA(HR,wtR,R,omega,A0,nq,wc,nw,kpt,alatt,blatt,nL)
% We build a supercell of nL unit cells, terminated vertical to b2 lattice
% Hence the first BZ is usually along x direction (a1 direction)
% This can be used to calculate surface state of a bullk system
% The orbital sequences are 1, ..., nw (in cell-1), 1, ..., nw (in cell-2)
% ..., until cell-nL.
% This nL unit cells form a principal layer for Green's function
nR = size(R,1);
R2max = max(R(:,2));
R2min = min(R(:,2));
nR2 = R2max-R2min+1;
HkqPn = zeros(nw,nw,nR2,nq);
Hk_qPn = zeros(nw,nw,nR2,nq);

FloPhase0 = FloPhase_ave(A0,R,wc,alatt,nw,nR,0);
Hk_q0Pn=buildHkAlongb1(R,HR.*FloPhase0,wtR,wc,nw,kpt,alatt,blatt);
% Hk_q0Pn is in size(nw,nw,nR2), same as dHdk_q0Pn (it is only dkx)
for iq = 1:nq
    FloPhaseq = FloPhase_ave(A0,R,wc,alatt,nw,nR,iq);
    FloPhase_q = FloPhase_ave(A0,R,wc,alatt,nw,nR,-iq);
    HkqPn(:,:,:,iq) = buildHkAlongb1(R,HR.*FloPhaseq,wtR,wc,nw,kpt,alatt,blatt);
    Hk_qPn(:,:,:,iq) = buildHkAlongb1(R,HR.*FloPhase_q,wtR,wc,nw,kpt,alatt,blatt);
end
HkPn = Hk_q0Pn;
for iq = 1:nq
    for iR2 = 1:nR2
    HkPn(:,:,iR2) = HkPn(:,:,iR2) ...
        + comm(Hk_qPn(:,:,iR2,iq),HkqPn(:,:,iR2,iq),-1)/omega/iq;
    end
end

H00 = zeros(nL*nw,nL*nw);
% V = zeros(nL*nw,nL*nw,3);
H01 = zeros(nL*nw,nL*nw);
R2max = max(R(:,2));
R2min = min(R(:,2));
for iL = 1:nL
for jL = 1:nL
    iR2 = jL-iL+1-R2min;
    iR2s = jL-iL+1-R2min-nL;
    if iR2 >=1 && iR2 <= R2max-R2min+1
        H00((iL-1)*nw+1:iL*nw,(jL-1)*nw+1:jL*nw) = HkPn(:,:,iR2);
        % V((iL-1)*nw+1:iL*nw,(jL-1)*nw+1:jL*nw,1) = dHdkPn(:,:,iR2,1);
    end
    if iR2s >= 1 && iR2s <= R2max-R2min+1
        H01((iL-1)*nw+1:iL*nw,(jL-1)*nw+1:jL*nw) = HkPn(:,:,iR2s);
    end
end
end

end
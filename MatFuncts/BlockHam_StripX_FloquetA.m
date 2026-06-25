
function [H00,H01] = BlockHam_StripX_FloquetA(HR,wtR,R,omega,A0,nq,wc,nw,kpt,alatt,blatt,nL)
% We build a supercell of nL unit cells, terminated vertical to a1
% The first BZ is usually along the b2 direction
% This can be used to calculate edge state of a 2D system
% The orbital sequences are 1, ..., nw (in cell-1), 1, ..., nw (in cell-2)
% ..., until cell-nL.
% This nL unit cells form a principal layer for Green's function
nR = size(R,1);
R1max = max(R(:,1));
R1min = min(R(:,1));
nR1 = R1max-R1min+1;
HkqPn = zeros(nw,nw,nR1,nq);
Hk_qPn = zeros(nw,nw,nR1,nq);
% dHdkqPn = zeros(nw,nw,nR1,nq);
% dHdk_qPn = zeros(nw,nw,nR1,nq);

FloPhase0 = FloPhase_ave(A0,R,wc,alatt,nw,nR,0);
% [Hk_q0Pn, dHk_q0Pn]=buildHkAlongb2(R,HR.*FloPhase0,wtR,wc,nw,kpt,alatt,blatt);
[Hk_q0Pn, ~]=buildHkAlongb2(R,HR.*FloPhase0,wtR,wc,nw,kpt,alatt,blatt);
% Hk_q0Pn is in size(nw,nw,nR2), same as dHdk_q0Pn (it is only dkx)
for iq = 1:nq
    FloPhaseq = FloPhase_ave(A0,R,wc,alatt,nw,nR,iq);
    FloPhase_q = FloPhase_ave(A0,R,wc,alatt,nw,nR,-iq);
    [HkqPn(:,:,:,iq), ~]...
        = buildHkAlongb2(R,HR.*FloPhaseq,wtR,wc,nw,kpt,alatt,blatt);
    [Hk_qPn(:,:,:,iq), ~]...
        = buildHkAlongb2(R,HR.*FloPhase_q,wtR,wc,nw,kpt,alatt,blatt);
end
HkPn = Hk_q0Pn;
% dHdkPn = dHk_q0Pn;
for iq = 1:nq
    for iR1 = 1:nR1
    HkPn(:,:,iR1) = HkPn(:,:,iR1) ...
        + comm(Hk_qPn(:,:,iR1,iq),HkqPn(:,:,iR1,iq),-1)/omega/iq;
    % dHdkPn(:,:,iR1) = dHdkPn(:,:,iR1) ...
    %     + comm(dHdk_qPn(:,:,iR1,iq),HkqPn(:,:,iR1,iq),-1)/omega/iq ...
    %     + comm(Hk_qPn(:,:,iR1,iq),dHdkqPn(:,:,iR1,iq),-1)/omega/iq;
    end
end

H00 = zeros(nL*nw,nL*nw);
% V = zeros(nL*nw,nL*nw);
H01 = zeros(nL*nw,nL*nw);
R1max = max(R(:,1));
R1min = min(R(:,1));
for iL = 1:nL
for jL = 1:nL
    iR1 = jL-iL+1-R1min;
    iR1s = jL-iL+1-R1min-nL;
    if iR1 >=1 && iR1 <= R1max-R1min+1
        H00((iL-1)*nw+1:iL*nw,(jL-1)*nw+1:jL*nw) = HkPn(:,:,iR1);
        % V((iL-1)*nw+1:iL*nw,(jL-1)*nw+1:jL*nw) = dHdkPn(:,:,iR1);
    end
    if iR1s >= 1 && iR1s <= R1max-R1min+1
        H01((iL-1)*nw+1:iL*nw,(jL-1)*nw+1:jL*nw) = HkPn(:,:,iR1s);
    end
end
end

end
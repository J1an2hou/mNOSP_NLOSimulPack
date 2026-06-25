function [H00,V,H01] = BlockHam_StripeY(HR,wtR,R,wc,nw,kpt,alatt,blatt,nL)
% We build a supercell of nL unit cells, terminated along b-lattice
% This can be used to calculate surface state of a bullk system
% The orbital sequences are 1, ..., nw (in cell-1), 1, ..., nw (in cell-2)
% ..., until cell-nL.
% This nL unit cells form a principal layer for Green's function
[HkPn,dHdkPn,~]=buildHkAlongb1(R,HR,wtR,wc,nw,kpt,alatt,blatt);
H00 = zeros(nL*nw,nL*nw);
V = zeros(nL*nw,nL*nw,3);
H01 = zeros(nL*nw,nL*nw);
R2max = max(R(:,2));
R2min = min(R(:,2));
for iL = 1:nL
for jL = 1:nL
    iR2 = jL-iL+1-R2min;
    iR2s = jL-iL+1-R2min-nL;
    if iR2 >=1 && iR2 <= R2max-R2min+1
        H00((iL-1)*nw+1:iL*nw,(jL-1)*nw+1:jL*nw) = HkPn(:,:,iR2);
        V((iL-1)*nw+1:iL*nw,(jL-1)*nw+1:jL*nw,1) = dHdkPn(:,:,iR2,1);
        V((iL-1)*nw+1:iL*nw,(jL-1)*nw+1:jL*nw,2) = dHdkPn(:,:,iR2,2);
        V((iL-1)*nw+1:iL*nw,(jL-1)*nw+1:jL*nw,3) = dHdkPn(:,:,iR2,3);
    end
    if iR2s >= 1 && iR2s <= R2max-R2min+1
        H01((iL-1)*nw+1:iL*nw,(jL-1)*nw+1:jL*nw) = HkPn(:,:,iR2s);
    end
end
end

end

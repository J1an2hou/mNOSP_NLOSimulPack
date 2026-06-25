function [H00,V,H01] = BlockHam_SurfZ(HR,wtR,R,wc,nw,kpt,alatt,blatt,nL)
% We build a supercell of nL unit cells, terminated along c-lattice
% This can be used to calculate surface state of a bullk system
% The orbital sequences are 1, ..., nw (in cell-1), 1, ..., nw (in cell-2)
% ..., until cell-nL.
% This nL unit cells form a principal layer for Green's function
[HkPn,dHdkPn,~]=buildHkAlongb1b2(R,HR,wtR,wc,nw,kpt,alatt,blatt);
H00 = zeros(nL*nw,nL*nw);
V = zeros(nL*nw,nL*nw,3);
H01 = zeros(nL*nw,nL*nw);
R3max = max(R(:,3));
R3min = min(R(:,3));
for iL = 1:nL
for jL = 1:nL
    iR3 = jL-iL+1-R3min;
    iR3s = jL-iL+1-R3min-nL;
    if iR3 >=1 && iR3 <= R3max-R3min+1
        H00((iL-1)*nw+1:iL*nw,(jL-1)*nw+1:jL*nw) = HkPn(:,:,iR3);
        V((iL-1)*nw+1:iL*nw,(jL-1)*nw+1:jL*nw,1) = dHdkPn(:,:,iR3,1);
        V((iL-1)*nw+1:iL*nw,(jL-1)*nw+1:jL*nw,2) = dHdkPn(:,:,iR3,2);
        V((iL-1)*nw+1:iL*nw,(jL-1)*nw+1:jL*nw,3) = dHdkPn(:,:,iR3,3);
    end
    if iR3s >= 1 && iR3s <= R3max-R3min+1
        H01((iL-1)*nw+1:iL*nw,(jL-1)*nw+1:jL*nw) = HkPn(:,:,iR3s);
    end
end
end

end

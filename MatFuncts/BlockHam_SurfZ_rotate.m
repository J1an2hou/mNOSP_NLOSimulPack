function [H00,V,H01] = BlockHam_SurfZ_rotate(U,HR,wtR,R,wc,nw,kpt,alatt,nL)
% We build a supercell of nL unit cells, terminated along c-lattice
% This can be used to calculate surface state of a bullk system
% The orbital sequences are 1, ..., nw (in cell-1), 1, ..., nw (in cell-2)
% ..., until cell-nL.
% This nL unit cells form a principal layer for Green's function
[HkPn,dHdkPn]=buildHkAlongb1b2_rotate(U,R,HR,wtR,wc,nw,kpt,alatt);
H00 = zeros(nL*nw,nL*nw);
V = zeros(nL*nw,nL*nw,3);
H01 = zeros(nL*nw,nL*nw);
Rallmax = max(max(R));
Rallmin = min(min(R));
for iL = 1:nL
for jL = 1:nL
    iR3 = jL-iL+1-Rallmin;
    iR3s = jL-iL+1-Rallmin-nL;
    if iR3 >=1 && iR3 <= Rallmax-Rallmin+1
        H00((iL-1)*nw+1:iL*nw,(jL-1)*nw+1:jL*nw) = HkPn(:,:,iR3);
        V((iL-1)*nw+1:iL*nw,(jL-1)*nw+1:jL*nw,1) = dHdkPn(:,:,iR3,1);
        V((iL-1)*nw+1:iL*nw,(jL-1)*nw+1:jL*nw,2) = dHdkPn(:,:,iR3,2);
        V((iL-1)*nw+1:iL*nw,(jL-1)*nw+1:jL*nw,3) = dHdkPn(:,:,iR3,3);
    end
    if iR3s >= 1 && iR3s <= Rallmax-Rallmin+1
        H01((iL-1)*nw+1:iL*nw,(jL-1)*nw+1:jL*nw) = HkPn(:,:,iR3s);
    end
end
end
if ~ishermitian(H00)
    warning('The Hamiltonian is not Hermitian')
end
H00 = roundn(H00,-6);
H01 = roundn(H01,-6);
V = roundn(V,-6);
end


function [Hk,dHdk] ...
        = buildHkAlongb1b2_rotate(U,R,HR,wtR,wc,nw,kpt,alatt)
% Here only R1 and R2 are Fourier transformed
% This is usually for surface spectral function calculation from 3D bulk
% Wannier functions

% This function perform a coordinate rotation of a given lattice and its
% fitted Wannier function. The U matrix is used to rotate the original
% lattice into a new set. The lattice and interacting R vector in Wannier
% function represetations are updated.
% According to the Wanniertools algorithm, the new corrdinate set is 
% x is along R1', z is along R1' x R2', y is along z \times y

% The U matrix has to be determined as input, which does not change the
% volume of the unit cell, and its components are 0 or 1 or -1
if (det(U)+1) < 1e-8
    U(3,:) = -U(3,:);
end
if abs(det(U) - 1) > 1e-8
    error('The rotation matrix U should not change the volume of unit cell')
end

dwc_rot = zeros(nw, 3);
nR = size(R,1);
Rallmax = max(max(R));
Rallmin = min(min(R));
nR3 = Rallmax-Rallmin+1;
Hk = zeros(nw,nw,nR3);
dHdk = zeros(nw,nw,nR3,3);
% d2Hdk2 = zeros(nw,nw,nR3,6);
dwc = zeros(nw,3);
for iw = 1:nw
    dwc(iw,:) = wc(iw,:) / alatt;
    dwc_rot(iw,:) = dwc(iw,:) / U;
end

for iR = 1:nR
    dwcdiff1 = dwc_rot(:,1)*ones(1,nw)-ones(nw,1)*dwc_rot(:,1)';
    dwcdiff2 = dwc_rot(:,2)*ones(1,nw)-ones(nw,1)*dwc_rot(:,2)';
    dwcdiff3 = dwc_rot(:,3)*ones(1,nw)-ones(nw,1)*dwc_rot(:,3)';
    R_rot = R(iR,:) / U;
    if norm(R_rot - round(R_rot)) < 1e-8
        R_rot = round(R_rot);
    else
        warning('Rotation of R lattice failed')
        disp(['R - rotated_R ', num2str(norm(R_rot - round(R_rot)))])
    end
    if R_rot(3) < Rallmin || R_rot(3) > Rallmax
        continue
    end
    aR_new = R * alatt / U;
    aR1 = aR_new(1);
    aR2 = aR_new(2);
    R1 = R_rot(1) - dwcdiff1;
    R2 = R_rot(2) - dwcdiff2;
    R3 = R_rot(3) - dwcdiff3;
    R3indx = R_rot(3)+1-Rallmin;
    dkdotdR = R1 * kpt(1) + R2*kpt(2) + R3*kpt(3);
    expkr = exp(1i*2*pi*dkdotdR);
    hr = reshape(HR(iR,:,:),nw,nw);
    Hk(:,:,R3indx) = Hk(:,:,R3indx) + wtR(iR) * hr .* expkr;
    dHdk(:,:,R3indx,1) = dHdk(:,:,R3indx,1)+1i*wtR(iR)*aR1.*hr.*expkr;
    dHdk(:,:,R3indx,2) = dHdk(:,:,R3indx,2)+1i*wtR(iR)*aR2.*hr.*expkr;
    % d2Hdk2(:,:,R3indx,1) = d2Hdk2(:,:,R3indx,1)-wtR(iR)*aR1.*aR1.*hr.*expkr;
    % d2Hdk2(:,:,R3indx,2) = d2Hdk2(:,:,R3indx,2)-wtR(iR)*aR2.*aR2.*hr.*expkr;
    % d2Hdk2(:,:,R3indx,3) = d2Hdk2(:,:,R3indx,6)-wtR(iR)*aR1.*aR2.*hr.*expkr;
end

rdn = -8;
Hk = roundn(Hk,rdn);
dHdk = roundn(dHdk,rdn);
% d2Hdk2 = roundn(d2Hdk2,rdn);
end
function [kptcell, nkpts] = k_mesh_BZ(nk1, nk2, nk3)
%GENERATE_KMESH Generate uniform k-mesh in first BZ centered at Gamma
%   b1, b2, b3: reciprocal lattice vectors (3x1 each)
%   nk1, nk2, nk3: number of k-points along each reciprocal direction
%
%   Output:
%       kcell: (Nk x 1) cell array, each cell contains [kx, ky, kz] in
%       direct coordinates (fraction of BZ vectors)
%
% Make 1D k vectors centered at Gamma
kx_vec = ((0:nk1-1) - floor(nk1/2))/nk1;
ky_vec = ((0:nk2-1) - floor(nk2/2))/nk2;
kz_vec = ((0:nk3-1) - floor(nk3/2))/nk3;

% Total number of k-points
nkpts = nk1 * nk2 * nk3;

% Preallocate cell
kptcell = cell(nkpts,1);

idx = 1;
for i1 = 1:nk1       % slowest index
    for i2 = 1:nk2
        for i3 = 1:nk3   % fastest index
            % Construct k-point in Cartesian coordinates
            kptcell{idx} =[kx_vec(i1), ky_vec(i2), kz_vec(i3)];
            idx = idx + 1;
        end
    end
end

%% Previous version
% Generate k coordinates in the first BZ
% kx = linspace(-0.5, 0.5, nk1+1);
% ky = linspace(-0.5, 0.5, nk2+1);
% kz = linspace(-0.5, 0.5, nk3+1);
% dkx = 1/nk1; dky = 1/nk2; dkz = 1/nk3;
% kx = kx+dkx/2; ky = ky+dky/2; kz = kz+dkz/2;
% [kx, ky, kz] = meshgrid(kx(1:nk1), ky(1:nk2), kz(1:nk3));
% kx = permute(kx,[2,1,3]);
% ky = permute(ky,[2,1,3]);
% kz = permute(kz,[2,1,3]);
% nkpts = nk1*nk2*nk3;
% kptcell = cell(nkpts,1);
% for ikpt = 1:nkpts
%     ikx = ceil(ikpt/(nk2*nk3)); % ikz is the fastest index
%     iky = ceil((ikpt-(ikx-1)*nk2*nk3)/nk3);
%     ikz = ceil(ikpt-(ikx-1)*nk2*nk3-(iky-1)*nk3);
%     kptcell{ikpt} = [kx(ikx,iky,ikz), ky(ikx,iky,ikz), kz(ikx,iky,ikz)];
% end

end
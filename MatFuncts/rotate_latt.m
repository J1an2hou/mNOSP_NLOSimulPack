function [aR_2dnew,R_new,wc_new,a_new,b_new] = rotate_latt(U,R,wc,alatt,blatt)
% This function perform a coordinate rotation of a given lattice and its
% fitted Wannier function. The U matrix is used to rotate the original
% lattice into a new set. The lattice and interacting R vector in Wannier
% function represetations are updated.
% According to the Wanniertools algorithm, the new corrdinate set is 
% x is along R1', z is along R1' x R2', y is along z \times y

% outputs:
% aR_2dnew is a rotated vector that is in the R1' x R2' plane in
% the new coordinate
% R_new is the direct R index (integers) that are in the new coordinate
% wc_new is the Wannier center absolute coordinates
% a_new and b_new are the real and reciprocal space vectors

% inputs:

% The U matrix has to be determined as input, which does not change the
% volume of the unit cell, and its components are 0 or 1 or -1
% R is the old direct R index (integers)
% wc is the Wannier centers in the old coordinate systems
% alatt and blatt are the real and reciprocal space vectors

if det(U) ~= 1
    error('The rotation matrix U should not change the volume of unit cell')
end
a_new = U * alatt;
b_new = blatt / U;
nR = size(R,1);
nw = size(wc,1);
wc_new = zeros(nw,3);
aR_2dnew = zeros(nR,3);
R_new = zeros(nR,3);
for iR = 1:nR
    R_idx = R(iR,:); % row vector of size (1,3), contains integers
    R_new(iR,:) = R_idx / U;
    r0 = R_new(iR,1) * a_new(1,:) + R_new(iR,2) * a_new(2,:); % this is a vector in the plane of R1' x R2'
    r1 = r0 / U;
    aR_2dnew(iR,:) = r1;
end

for iw = 1:nw
    wc_new(iw,:) = wc(iw,:) / U;
end

end
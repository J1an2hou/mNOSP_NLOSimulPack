function [alatt, blatt, dxyz, mxyz, nat] = TB_readlattice(filename)
% This file is quite similar as POSCAR in VASP
% Line 0: a file head
% Line 1: scaling factor
% Line 2: R1
% Line 3: R2
% Line 4: R3
% Line 5: number of sites
% Line 6: a dummy line
% Line 7: coordinates, relative to lattice
% ......
% Line xx: spin polarization of site 1, Cartesian coordinate
% Line xx+1: spin polarization of site 2, Cartesian coordinate
% ......
fileID = fopen(filename);
fgetl(fileID);
scale = sscanf( fgetl(fileID), '%f' );
alatt = zeros(3,3);
for i = 1:3
  alatt(i,:) = sscanf( fgetl(fileID), '%f' );
end
alatt = scale * alatt;
V = det(alatt);
fgetl(fileID);
nat = sscanf( fgetl(fileID), '%f');
nat = sum(nat);
fgetl(fileID);
dxyz = zeros(nat,3);
mxyz = zeros(nat,3);
for i = 1:nat
    dxyz(i,:) = sscanf( fgetl(fileID), '%f' );
end
fgetl(fileID); % A line separating coordinates and spin polarization
for i = 1:nat
    mxyz(i,:) = sscanf(fgetl(fileID),'%f');
end
fclose(fileID);
blatt = zeros(3,3);
blatt(1,:) = 2*pi*cross(alatt(2,:),alatt(3,:))/V;
blatt(2,:) = 2*pi*cross(alatt(3,:),alatt(1,:))/V;
blatt(3,:) = 2*pi*cross(alatt(1,:),alatt(2,:))/V;

end


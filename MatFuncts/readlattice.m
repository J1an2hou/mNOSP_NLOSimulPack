function [alatt, blatt, Vlatt] = readlattice(filename)
% Reading lattice constants and evaluating reciprocal lattice
% Format copies from VASP-POSCAR
file = fopen(filename);
fgetl(file); 
scale = sscanf( fgetl(file), '%f' );
alatt = zeros(3,3);
for i = 1:3
  alatt(i,:) = sscanf( fgetl(file), '%f' );
end
alatt = scale * alatt;
Vlatt = det(alatt);

blatt = zeros(3,3);
blatt(1,:) = 2*pi*cross(alatt(2,:),alatt(3,:))/Vlatt;
blatt(2,:) = 2*pi*cross(alatt(3,:),alatt(1,:))/Vlatt;
blatt(3,:) = 2*pi*cross(alatt(1,:),alatt(2,:))/Vlatt;
fclose(file);
end
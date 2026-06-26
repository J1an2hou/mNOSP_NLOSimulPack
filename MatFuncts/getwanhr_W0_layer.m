
function [R,HR,wtR,wc,wU0,wl,nw,nR] = getwanhr_W0_layer(seedname)
file = fopen(strcat([seedname, '_hr.dat']), 'r');
fgetl(file);
line = fgetl(file);
nw = sscanf(line, '%d');
line = fgetl(file);
nR = sscanf(line, '%d');

Rdeg = [];
Rlines = ceil(nR/15); % per the format of wannier90_hr, 15 entries per line
for iline = 1:Rlines
	line = fgetl(file);
  Rdeg = [Rdeg, sscanf(line, '%d')'];
end
wtR = 1./Rdeg;

C = textscan(file,'%d %d %d %d %d %n %n');
R1 = C{1};
R2 = C{2};
R3 = C{3};
Hr1 = C{6};
Hr2 = C{7};

fclose(file);
R = zeros(nR,3);
HR = zeros(nR,nw,nw);
for i = 1:nR
    R(i,1) = R1((i-1)*nw^2+1);
    R(i,2) = R2((i-1)*nw^2+1);
    R(i,3) = R3((i-1)*nw^2+1);
    for j = 1:nw
        for k = 1:nw
            HR(i,j,k) = Hr1((i-1)*nw^2+(k-1)*nw+j) ...
                +1i*Hr2((i-1)*nw^2+(k-1)*nw+j);
        end
    end
    
end
% symmetrize HR
for i = 1:nR
  HR(i,:,:) = 0.5*( HR(i,:,:) + conj(permute(HR(nR+1-i,:,:), [1,3,2])) );
end
% wc = zeros(nw,3);

fileID = fopen(strcat([seedname,'_centres-W0-layer.xyz']),'r');
fgetl(fileID);
fgetl(fileID);
C = textscan(fileID,'%*s %f %f %f %f %f');
wc1 = C{1};
wc2 = C{2};
wc3 = C{3};
U01 = C{4}; % if giving a disorder or not
proj_layer = C{5};
fclose(fileID);
wc = zeros(nw,3);
wU0 = zeros(nw,1);
wl = zeros(nw, 1);
for i = 1:nw
    wc(i,1) = wc1(i);
    wc(i,2) = wc2(i);
    wc(i,3) = wc3(i);
    wU0(i) = U01(i);
    wl(i) = proj_layer(i);
end

end


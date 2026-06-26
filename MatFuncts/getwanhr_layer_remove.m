function [R,HR,wtR,wc,wl,xl,nw,nR] = getwanhr_layer_remove(seedname)
% Read the Wannier (seedname usually wannier90) results
% Input files include wannier90_hr.dat, (modified) wannier90_centres.xyz
filehr = fopen(strcat([seedname, '_hr.dat']), 'r');
fgetl(filehr);
line = fgetl(filehr);
nw = sscanf(line, '%d');
line = fgetl(filehr);
nR = sscanf(line, '%d');

Rdeg = [];
Rlines = ceil(nR/15); % per the format of wannier90_hr, 15 entries per line
for iline = 1:Rlines
	line = fgetl(filehr);
  Rdeg = [Rdeg, sscanf(line, '%d')'];
end
wtR = 1./Rdeg;

C = textscan(filehr,'%d %d %d %d %d %n %n');
R1 = C{1};
R2 = C{2};
R3 = C{3};
Hr1 = C{6};
Hr2 = C{7};

fclose(filehr);
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

filexyz = fopen(strcat([seedname,'_centres-layers-remove.xyz']),'r');
fgetl(filexyz);
fgetl(filexyz);
C = textscan(filexyz,'%s %f %f %f %f');
Del = C{1};
wc1 = C{2};
wc2 = C{3};
wc3 = C{4};
layer = C{5};
fclose(filexyz);
wc = zeros(nw,3);
wl = zeros(nw,1);
xl = zeros(nw,1);
for i = 1:nw
    if string(Del(i)) == 'N'
        xl(i) = -1;
    end
    wc(i,1) = wc1(i);
    wc(i,2) = wc2(i);
    wc(i,3) = wc3(i);
    wl(i) = layer(i);
end

end

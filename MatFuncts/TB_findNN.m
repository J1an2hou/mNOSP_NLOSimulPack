function [delta, jindex] = TB_findNN(alatt,nat,nc,d,dxyz)
dthres = 0.1; % bond length measurement threshold
xyz = zeros(nat,3);
for i = 1:nat
    xyz(i,:) = dxyz(i,:) * alatt;
end
delta = zeros(nat,nc,3); % distancing vectors from atom-i to j
jindex = zeros(nat,nc);
for i = 1:nat % home site
    nnc = 0;
for j = 1:nat
    for nR1 = -1:1
    for nR2 = -1:1
    for nR3 = -1:1
        drij = xyz(j,:)  - xyz(i,:)...
            + nR1*alatt(1,:)+nR2*alatt(2,:)+nR3*alatt(3,:);
        dist = sqrt(drij * drij');
        if abs(dist - d) < dthres
            nnc = nnc+1;
            delta(i,nnc,:) = drij;
            jindex(i,nnc) = j;
        end
    end
    end
    end
end
    if nnc > nc
        disp(['Calculated coordinate ',num2str(nnc),', while given coordinate number ',num2str(nc)])
        error('Check the coordinate number')
    end
end
end
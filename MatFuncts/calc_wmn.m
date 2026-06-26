function w = calc_wmn(wavefunct, d2Hdk2,nw)

w = zeros(nw,nw,3,3);

for i = 1:3
    for j = 1:3
        k = voigt(i,j);
        w(:,:,i,j) = wavefunct' * d2Hdk2(:,:,k) * wavefunct;
    end
end
end


function lpmn = calc_proj(wavefunct, lproj)
nlayer = size(lproj,3);
nw = size(lproj,1);
lpmn = zeros(nw, nw, nlayer);
for ilayer = 1:nlayer
    lpmn(:,:,ilayer) = wavefunct' * lproj(:,:,ilayer) * wavefunct;
end
end
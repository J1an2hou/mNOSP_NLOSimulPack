function vmn = calc_vmn(wavefunct, dHdk, nw)
vmn = zeros(nw, nw, 3);

vmn(:,:,1) = wavefunct' * dHdk(:,:,1) * wavefunct;
vmn(:,:,2) = wavefunct' * dHdk(:,:,2) * wavefunct;
vmn(:,:,3) = wavefunct' * dHdk(:,:,3) * wavefunct;
end
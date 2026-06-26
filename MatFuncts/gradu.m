function dudk = gradu(wavefunct, vmn, omn, nw, ksi)
dudk = zeros(nw, nw, 3);
for dir = 1:3
    vovero = - vmn(:,:,dir) .* omn ./ (omn.^2 + ksi^2);
    dudk(:,:,dir) = wavefunct * vovero;
end
end
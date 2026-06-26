function [xmn, ymn, zmn] = calc_position(wavefunct, wc)
x = wc(:,1);
y = wc(:,2);
z = wc(:,3);
x_ave = mean(x);
y_ave = mean(y);
z_ave = mean(z);
x = diag(x - x_ave);
y = diag(y - y_ave);
z = diag(z - z_ave);
xmn = wavefunct' * x * wavefunct;
ymn = wavefunct' * y * wavefunct;
zmn = wavefunct' * z * wavefunct;
end
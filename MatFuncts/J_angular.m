function [Jx, Jy, Jz] = J_angular(j)
% transform the angular momentum matrix in spherical harmonics
% basis to | j, j-1, ..., -j > basis
 
% we define U by | atomic orbitals > = U | j, j-1, ..., -j >
% therefore, J_|atomicorbitals> = U * J_sp * inv(U);
if j == 0      % | s >
  U = eye(1);
elseif j == 1  % | pz, px, py >
  U = sqrt(1/2) * [0,  sqrt(2), 0;
                   -1,   0,     1;
                   1i,   0,    1i];
elseif j == 2  % | dz2, dxz, dyz, dx2-y2, dxy >
  U = sqrt(1/2) * [0,  0,  sqrt(2),  0,  0;
                   0, -1,    0,      1,  0;
                   0, 1i,    0,     1i,  0;
                   1,  0,    0,      0,  1;
                  -1i, 0,    0,      0, 1i];
else
  error('Angular Momentum should be 0, 1 or 2');
end
[Jx_sp, Jy_sp, Jz_sp] = J_SphericalHarmonics(j);
Jx = U * Jx_sp / U;
Jy = U * Jy_sp / U;
Jz = U * Jz_sp / U;
Jx = roundn(Jx,-10);
Jy = roundn(Jy,-10);
Jz = roundn(Jz,-10);
% We see a sign problem (for Jx and Jz) as compared with previous works
% PRB 98, 214405 (2018); PRB 103, 085113 (2011); PRB 77, 165117 (2008)
% so we add a sign
Jx = -Jx; Jz = -Jz;
end


% We incorporate linear and nonlinear optical responses here
% All calculations are based on Wannier functions
% -----------------------------------------------------------------------
% The following papers need to be cited if you use our code
% (dielectrics) J. Zhou, H. Xu, Y. Shi, J. Li, Adv. Sci. 2021, 8, 2003832
% (Hall effect) Y. Pan, J. Zhou, Phys. Rev. Appl. 2020, 14, 014024
% (BPV effect) X. Mu, Y. Pan, J. Zhou, npj Comput. Mater. 2021, 7, 61
% (Edelstein effect) Y. Sun, X. Mu, Q. Xue, J. Zhou, Adv. Opt. Mater. 2022, 10, 2200428
% (Floquet engineering) to be announced
% -----------------------------------------------------------------------
% Author: Jian Zhou
%% Starting
clear;
ncores = 56; % number of cores for parallel processing (or use parpool('Threads') for automatic pool)
%% Physical Parameters
pie3h = 765.1;
% This is constant conversion from pi*e^3/hbar to \mu / V^2
pie3h2 = 1.16; % To convert pi*e^3/hbar^2 to A/V^2/ps

%% Read Files
[R,HR,wtR,wc,wl,nw,nR] = getwanhr_layer('wannier90');
% [R,HR,wtR,wc,nw,nR] = getwanhr('wannier90');
[alatt, blatt, Vlatt] = readlattice('POSCAR');
function wk = spectral_weight_unfolding(wave_sc, wc_sc, nw_sc, bG_sc, N_sc)
% This script evaluates the spectral weight of each state at a specific k
% point, according to Eq. (35) of EPL, 107 (2014) 27006.
% We assume that the basis set of a supercell is duplicated by unit cell
% basis set, and multiplied by N_sc times. In this case, one can figure out
% the equivalent orbitals in the supercell.
% In the equation, one needs to evaluate the \sum_j^N c^{a+j}(k)
% In this case, we first multiply phase factor e^{i*G*R_{a+j}}, which is
% denoted as expkr in the following. Then, we convolute this summation into
% a new wave function, that is in the size of (nw_uc, uw_sc), using a
% conv_M matrix
%
% Parameter sizes:
% bG_sc is in size (1 x 3)
% wave_sc is in size of (nw_sc * nw_sc), with each column gives
% wavefunctions
% Output: wk is in size of (nw_sc * 1), giving spectral weight of each k

nw_uc = nw_sc / N_sc;
expkr = exp(-1i* wc_sc * bG_sc'); % size of (nw_sc x 1)
wave_pbz = expkr .* wave_sc;
% Build a converting matrix
conv_M = kron(ones(1, N_sc), eye(nw_uc, nw_uc));
wave_sum = conv_M * wave_pbz;
wk = real(diag(wave_sum' * wave_sum));
wk = wk/N_sc;
end
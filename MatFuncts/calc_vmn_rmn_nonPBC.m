function [vmn_i,rmn_i] = calc_vmn_rmn_nonPBC(omn,wvf,wc,direct)
% We evaluate the < m | r | n > = wvf' * wc * wvf for position operator
% Then the velocity operator is calculated according to
% vmn = i * omn .* rmn
% Both vmn_i and rmn_i are in size of (nw x nw)
% Note that here the vmn and rmn are in the direction of i (direct)
r = wc(:,direct);
r_ave = mean(r);
r = r - r_ave;
rmn_i = wvf' * diag(r) * wvf;
vmn_i = 1i*omn .* rmn_i;
end
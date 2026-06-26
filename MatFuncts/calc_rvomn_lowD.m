function [r,v,o] = calc_rvomn_lowD(v,energy,nw,nonperiod_dir,wvf,wc,eta)
% We calculate the low dimensional vmn and rmn with some non-periodic
% boundary, the nonperiod_dir indicates which direction (1 - 3) should be
% replaced by the nonperiodic position and velocity
% The nonperiod_dir must be a list, or a number showing the nonperiodic
% direction. In a lot of cases, it is 3 (z direction for 2D material)
[xmn, ymn, zmn] = calc_position(wvf, wc);
r_pos = zeros(nw,nw,3);
r_pos(:,:,1) = xmn;
r_pos(:,:,2) = ymn;
r_pos(:,:,3) = zmn;
r = zeros(nw,nw,3);
o = energy*ones(1,nw) - ones(nw,1)*energy';
for idir = 1:3
    if ismember(idir, nonperiod_dir) % nonperiodic direction
        r(:,:,idir) = r_pos(:,:,idir);
        v(:,:,idir) = 1i*r(:,:,idir) .* o; % update velocity operator
    else
        r(:,:,idir) = -1i*v(:,:,idir) .* o ./ (o.^2+eta^2);
    end
end

% v1 = v(:,:,1);
% v2 = v(:,:,2);
% v3 = v(:,:,3);
% r1 = -1i*v1 .* o ./ (o.^2 + eta^2);
% r2 = -1i*v2 .* o ./ (o.^2 + eta^2);
% r3 = -1i*v3 .* o ./ (o.^2 + eta^2);
% r(:,:,1) = r1;
% r(:,:,2) = r2;
% r(:,:,3) = r3;
end
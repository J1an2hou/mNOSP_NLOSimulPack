function nu = TB_nu_chirality(nat, nc1, delta1, jindex1, nc2, delta2, jindex2)
% compute all the chirality for soc term
% The nc1, delta1, and jindex1 are the parameters for NN bonds
% The nc2, delta2, and jindex2 can be either NN or NNN, depending on which
% type of soc chirality you are working on
if nc2 ~= size(delta2,2) || nc1 ~= size(delta1,2)
    error('Wrong match between NN or NNN size parameters')
end
nu = zeros(nat, nc2);
for iat = 1:nat
    for ic2 = 1:nc2
        jat = jindex2(iat,ic2);
        if jat == 0
            continue
        end
        delta_ij = reshape(delta2(iat,ic2,:),1,3); % from i to j (over k)
        nu(iat, ic2) = nu_ikj(iat, jat, nc1, delta1, jindex1, delta_ij);
    end
end

end

function nu_ij = nu_ikj(iat, jat, nc, delta, jindex, delta_ij)
% ==========================================================
% Explicitly compute nu_ij by finding the intermediate atom k
% for intrinsic kagome SOC
% It can be used for either NN or NNN bonding, when one can find a k
% to connect i and j atoms
% ==========================================================
tol = 1e-4;
if nc ~= size(delta,2)
    error('Wrong match between delta and nc parameter')
end
% loop over NN k of i
for ic1 = 1:nc
    kat = jindex(iat,ic1);
    if kat == 0
        continue
    end
    dik = reshape(delta(iat,ic1,:),1,3);
    % check if k connects to j by NN
    for jc1 = 1:nc
        if jindex(kat,jc1) ~= jat
            continue
        end
        dkj = reshape(delta(kat,jc1,:),1,3);
        dijki = dik + dkj - delta_ij;
        dist = sqrt(sum(dijki.^2));
        % disp(dist)
        if dist < tol
            % chirality from cross product
            chir = dik(1)*dkj(2) - dik(2)*dkj(1);
            nu_ij = sign(chir);
            return
        end
    end
end

error('No valid i -> k -> j path found');
end
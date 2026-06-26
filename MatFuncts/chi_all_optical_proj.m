function [chiee_k,chiem_k,chime_k,chimm_k] ...
    = chi_all_optical_proj(o,f,r,p,m,hv,ksi)
% This function evaluates the Kubo formula for linear optical effect
% chiee - traditional optical dielectric function
% chiem - optical magnetoelectric coupling
% chime - inverse optical ME coupling
% chimm - magnetic "dielectric" function, usually does not take part in response
% In the input matrix, they are
% o     - energy difference matrix
% f     - occupation difference matrix
% r     - position matrix
% p     - projection matrix
% m     - spin / orbital matrix, with unit of mu_B (if spin, maybe g-factor of 2 is needed)
% hv    - incident photon energy
% ksi   - broadening factor
nhv = length(hv);
nwl = size(p,3);
chiee_k = zeros(3, 3, nhv,nwl);
chiem_k = zeros(3, 3, nhv,nwl);
chime_k = zeros(3, 3, nhv,nwl);
chimm_k = zeros(3, 3, nhv,nwl);

for iwl = 1:nwl
    proj = p(:,:,iwl);
    for ihv = 1:nhv
        omega = hv(ihv);
        omn_omega = o-omega-1i*ksi;
        for ii = 1:3
        for jj = 1:3
        % chiee_k(ii,jj,ihv,iwl) = trace(((r(:,:,ii)*proj) ./ omn_omega) * (f .* r(:,:,jj)));
        % chiem_k(ii,jj,ihv,iwl) = trace(((r(:,:,ii)*proj) ./ omn_omega) * (f .* m(:,:,jj)));
        % chime_k(ii,jj,ihv,iwl) = trace(((m(:,:,ii)*proj) ./ omn_omega) * (f .* r(:,:,jj)));
        % chimm_k(ii,jj,ihv,iwl) = trace(((m(:,:,ii)*proj) ./ omn_omega) * (f .* m(:,:,jj)));
        chiee_k(ii,jj,ihv,iwl) = trace(((comm(r(:,:,ii),proj,1)/2) ./ omn_omega) * (f .* r(:,:,jj)));
        chiem_k(ii,jj,ihv,iwl) = trace(((comm(r(:,:,ii),proj,1)/2) ./ omn_omega) * (f .* m(:,:,jj)));
        chime_k(ii,jj,ihv,iwl) = trace(((comm(m(:,:,ii),proj,1)/2) ./ omn_omega) * (f .* r(:,:,jj)));
        chimm_k(ii,jj,ihv,iwl) = trace(((comm(m(:,:,ii),proj,1)/2) ./ omn_omega) * (f .* m(:,:,jj)));
        end
        end
    end
end


end

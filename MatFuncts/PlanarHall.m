function [sigma_sp,sigma_orb,bc_layer] ...
    = PlanarHall(v,r,o,en,p,sp1,sp2,orb1,orb2,B_phi,cp,EF,ksi)
% This function is to evaluate the first-order corection of planar Hall
% effect under an in-plane B field (1 Tesla). The zero-order is trivial anomalous
% Hall conductance that comes from (intrinsically) Berry curvature, while
% the first-order correction is expressed as
% \sigma^(1) = \int_k (df/dE) * [v_n \times A_n + m_n * Omega_n]\cdot B
% Here, intraband A_n = -2*\hbar*Im(v_nl * m_ln)/(e_n - e_l}^3
% For simplicity, only sigma_xy component is evaluated, separately for spin
% moment and orbital moment responses
% The equations are taking from Nano Lett. 25, 10096 (2025), Eq. (2) and (S10)
ncp = length(cp);
nwl = size(p,3);
sigma_sp = zeros(ncp, nwl);
sigma_orb = zeros(ncp, nwl);
bc_layer = zeros(ncp, nwl);

Bx = cos(B_phi); By = sin(B_phi);
m_orb = orb1 * Bx + orb2 * By;
m_sp = sp1 * Bx + sp2 * By;
v1 = v(:,:,1);
v2 = v(:,:,2);
r1 = r(:,:,1);
r2 = r(:,:,2);

mt_orb = m_orb .* o ./ (o.^2+ksi^2);
mt_sp = m_sp .* o ./ (o.^2+ksi^2);
mn_orb = real(diag(m_orb));
mn_sp = real(diag(m_sp));

for icp = 1:ncp
    mu = EF + cp(icp);
    dfde = -delta_funct(en - mu, ksi);
    f = 1./(exp(en - mu)/ksi + 1);
    for iwl = 1:nwl
    pj = p(:,:,iwl);
    v1_proj = comm(v1, pj, 0.5);
    v1np = diag(real(v1_proj));
    v2n = diag(real(v2));
    r1_proj = comm(r1, pj, 0.5);
    Omega_n = diag(imag(r1_proj*r2 - r2*r1_proj));
    sigma_sp(icp, iwl) = sum(2*(v1np .* real(diag(r2 * mt_sp)) - v2n .* real(diag(r1_proj * mt_sp))) .* dfde ...
        + (mn_sp.*Omega_n).* dfde);
    sigma_orb(icp, iwl) = sum(2*(v1np .* real(diag(r2 * mt_orb)) - v2n .* real(diag(r1_proj * mt_orb))) .* dfde ...
        + (mn_orb.*Omega_n).* dfde);
    bc_layer(icp, iwl) = sum(f .* Omega_n);
    end

end
sigma_sp = sigma_sp * 2; % spin Lande g-factor is 2
sigma_sp = sigma_sp *   5.788e-5; % muB * Tesla to eV
sigma_orb = sigma_orb * 5.788e-5;
end
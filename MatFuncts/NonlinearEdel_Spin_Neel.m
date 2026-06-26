function [sp_even,sp_odd,nl_even,nl_odd,gkk,gsk,gnk] ...
    = NonlinearEdel_Spin_Neel(v,o,en,sp1,sp2,sp3,nl1,nl2,nl3,nw,cp,EF,direct,ksi)
% We evaluate the time-even and time-odd second order Edelstein effect
% responses, according to PRL 2022, 129, 086602 and PRL 2023, 130, 166302
% In order to do them, we need to evaluate the Berry connection
% polarizability and the (spin-k) mixed Berry connection polarizability
% The BCP denoted as G_kk and (spin-k) mixed BCP denoted as G_sk
% The anomalous spin correction is denoted as R_sk
% The equation is e^2*hbar^2 * int_{2D k} dk * v * v * S / (\hbar*omn)^4
% The integration in 2D k space will give the unit metrices with muB/V^2

% e_charge = -1;
ncp = length(cp);
ndir = size(direct,1);
sp_even = zeros(ndir, ncp);
sp_odd = zeros(ndir, ncp);
nl_even = zeros(ndir, ncp);
nl_odd = zeros(ndir, ncp);
gkk = zeros(nw, ndir); % berry connection polarizability, BCP
gsk = zeros(nw, ndir); % mixed BCP (spin and k)
gnk = zeros(nw, ndir); % mixed BCP (neel and k)

for idir = 1:ndir
    a = direct(idir, 1);
    b = direct(idir, 2);
    c = direct(idir, 3);
    va = v(:,:,a);
    vb = v(:,:,b);
    va_diag = real(diag(va));
    vb_diag = real(diag(vb));
    if c == 1
        sp = sp1; nl = nl1;
    elseif c == 2
        sp = sp2; nl = nl2;
    elseif c == 3
        sp = sp3; nl = nl3;
    else
        error('Wrong spin or Neel vector direction')
    end
    o2 = o.^2;
    o3 = o.^3;
    o4 = o.^4;

    % Fermi surface contributions
    G_kkab = 2 * diag((va .* o3 ./ (o3.^2+ksi^3)) * vb);
    % G_kkba = 2 * diag((vb .* o3 ./ (o3.^2+ksi^3)) * va);
    G_skca = -2 * diag((sp .* o3 ./ (o3.^2+ksi^3)) * va);
    G_skcb = -2 * diag((sp .* o3 ./ (o3.^2+ksi^3)) * vb);
    G_nkca = -2 * diag((nl .* o3 ./ (o3.^2+ksi^3)) * va);
    G_nkcb = -2 * diag((nl .* o3 ./ (o3.^2+ksi^3)) * vb);
    R_skca = 2 * diag((sp .* o2 ./ (o2.^2+ksi^2)) * va);
    R_skcb = 2 * diag((sp .* o2 ./ (o2.^2+ksi^2)) * vb);
    R_nkca = 2 * diag((nl .* o2 ./ (o2.^2+ksi^2)) * va);
    R_nkcb = 2 * diag((nl .* o2 ./ (o2.^2+ksi^2)) * vb);
    G_kkab = real(G_kkab);
    % G_kkba = real(G_kkba);
    G_skca = real(G_skca);
    G_skcb = real(G_skcb);
    G_nkca = real(G_nkca);
    G_nkcb = real(G_nkcb);
    gkk(:, idir) = G_kkab;
    gsk(:, idir) = G_skca;
    gnk(:, idir) = G_nkca;
    R_skca = imag(R_skca);
    R_skcb = imag(R_skcb);
    R_nkca = imag(R_nkca);
    R_nkcb = imag(R_nkcb);
    sc = real(diag(sp));
    nc = real(diag(nl));

    % Fermi sea contributions
    sc_diff = sc * ones(1,nw) - ones(nw,1) * sc';
    nc_diff = nc * ones(1,nw) - ones(nw,1) * nc';
    s_on_o = sp .* o ./ (o.^2 + ksi);
    n_on_o = nl .* o ./ (o.^2 + ksi);
    va_on_o3 = va .* o3 ./ (o3.^2+ksi^3);
    vb_on_o3 = vb .* o3 ./ (o3.^2+ksi^3);
    facA_sp = -3* sc_diff .* va.' .* vb .* o4 ./ (o4.^2 + ksi^4);
    facA_nl = -3* nc_diff .* va.' .* vb .* o4 ./ (o4.^2 + ksi^4);
    facB_sp = vb_on_o3 .* (s_on_o * va).' + va_on_o3 .* (s_on_o * vb).';
    facB_nl = vb_on_o3 .* (n_on_o * va).' + va_on_o3 .* (n_on_o * vb).';
    facC_sp = (s_on_o * va) .* vb_on_o3.' + (s_on_o * vb) .* va_on_o3.';
    facC_nl = (n_on_o * va) .* vb_on_o3.' + (n_on_o * vb) .* va_on_o3.';
    facA_sp = real(facA_sp);
    facA_nl = real(facA_nl);
    facB_sp = real(facB_sp);
    facB_nl = real(facB_nl);
    facC_sp = real(facC_sp);
    facC_nl = real(facC_nl);
    % End of Fermi sea contributions

    for icp = 1:ncp
    mu = EF + cp(icp);
    f = 1./(1+exp((en-mu)/ksi));
    f_diff = f*ones(1,nw) - ones(nw,1) * f';
    dfde = delta_funct(en - mu, ksi);
    
    sp_odd(idir, icp) = -sum((sc .* G_kkab + vb_diag .* G_skca + va_diag .* G_skcb) .* dfde)/2 ...
        + trace(f_diff * (facA_sp + facB_sp - facC_sp))/2;
    sp_even(idir, icp) = -sum((va_diag .* R_skcb + vb_diag .* R_skca) .* dfde);
    nl_odd(idir, icp) = -sum((nc .* G_kkab + vb_diag .* G_nkca + va_diag .* G_nkcb) .* dfde)/2 ...
        + trace(f_diff * (facA_nl + facB_nl - facC_nl))/2;
    nl_even(idir, icp) = -sum((va_diag .* R_nkcb + vb_diag .* R_nkca) .* dfde);
    end

end
% No converting units. The time-even is in unit of mu_B / V^2 /ps, and the
% time-odd will be in  unit of mu_B / V^2 /ps
end


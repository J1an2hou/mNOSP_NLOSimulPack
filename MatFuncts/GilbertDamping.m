function alpha_on_site_sf = GilbertDamping(H,p,EF,Enk,nw,wf,m_dir,gamma, T)
% This script calculates the intrinsic Gilbert damping alpha
% Then the alpha term is evaluated for different sublattices separately
% The formula is alpha = \sum_{m,n} < n | [sigma^A, Hsoc] | m > * < m
% |[sigma^B, H_soc] | n > * A_n(EF) * A_m(EF)
% Here A_n(EF) is the spectral function of band-n referenced to Fermi
% energy
% Inputs are here
% H       - Hamiltonian that is a matrix of nw x nw
% p        - usual projection operator, not multiplied with wavefunct
% EF      - Fermi energy
% Enk    - eigenenergies for a k point, in size of (nw, 1)
% nw     - size of Hamiltonian
% wf      - wavefunction, in size of (nw, nw)
% m_dir - direction of the equilibrium magnetization, size of (nwl, 3)
% gamma - Damping factor for Lorentz expansion of delta function
% T         - Temperature (optional, if not indicated, determined by gamma)
kB = 8.617333262145e-5;  % eV/K

nksi = length(gamma);
Lande_g = 2; % Lande g-factor
nwl = size(p, 3);
sx = Pauli_mat(nw, 1);
sy = Pauli_mat(nw, 2);
sz = Pauli_mat(nw, 3);
m_atm = size(m_dir, 1);
if m_atm == 1 && nwl > 1
    m_dir = kron(ones(nwl, 1), m_dir);
end
% alpha_tensor_sf = zeros(3, 3, nwl,nksi); % Fermi surface
alpha_on_site_sf = zeros(nwl,nksi); % Fermi surface
for iksi = 1:nksi
    ksi_loop = gamma(iksi);
    for iwl = 1:nwl
        pa = p(:,:,iwl);
        [az, el, ~] = cart2sph(m_dir(iwl,1), m_dir(iwl,2), m_dir(iwl,3));
        ma_the = pi/2 - el;
        ma_phi = mod(az, 2*pi);
        s_perp1 = sx*cos(ma_the)*cos(ma_phi) + sy*cos(ma_the)*sin(ma_phi) - sz*sin(ma_the);
        s_perp2 = -sx * sin(ma_phi) + sy * cos(ma_phi);
        stor_perp1_proj = comm(comm(s_perp1, H,-1),pa,1)/2;
        stor_perp2_proj = comm(comm(s_perp2, H,-1),pa,1)/2;
        stor_perp1 = comm(s_perp1, H,-1);
        stor_perp2 = comm(s_perp2, H,-1);
        stor1pmn = wf' * stor_perp1_proj * wf;
        stor2pmn = wf' * stor_perp2_proj * wf;
        stor1mn = wf' * stor_perp1 * wf;
        stor2mn = wf' * stor_perp2 * wf;
        % Wmn = inter_broadening(Enk,EF,nw,ksi_loop);
            if nargin == 8 % temperature not specified
                T = ksi_loop / kB;
            end
        Inte_fac = Two_Spectral_Overlap(Enk, ksi_loop, EF, T);
        alpha_on_site_sf(iwl,iksi) = -trace( Inte_fac ...
            * (stor1pmn .* stor1mn.' + stor2pmn .* stor2mn.'));
        % alpha_on_site_sf(iwl,iksi) = trace( Inte_fac ...
        %     * (abs(stor1pmn).^2 + abs(stor2pmn).^2));
    end
end

alpha_on_site_sf = pi * Lande_g * real(alpha_on_site_sf);

end


function W = inter_broadening(En,EF,nw,ksi)
% This function evaluates the W_nm matrix for interband transition
% This W_nm = pi * (fn-fm) / (En-Em)^2 but with a broadening factor ksi

atan_En = atan((En - EF)/ksi);
atan_diffE = atan_En * ones(1, nw) - ones(nw, 1)*atan_En';
diffE = En * ones(1,nw)-ones(nw,1)*En';
Abar = (En - EF) ./ ((En - EF).^2 + ksi^2);
A = pi*delta_funct(En-EF,ksi);
W = -atan_diffE .* (diffE.^2) ./ (diffE.^4 + ksi^4) ...
    +ksi * (A * Abar' + Abar * A') .* diffE ./ (diffE.^2 + ksi^2);

% alternatively in a simple way
% fn = 1./(exp((En-EF)/ksi) + 1);
% fdiff = fn * ones(1,nw)-ones(nw,1)*fn';
% W = pi * fdiff .* (diffE.^2) ./ (diffE.^4 + ksi^4);
end

function Imn = Two_Spectral_Overlap(Enk, eta, EF, T)
%  
%
% Inputs:
%   Enk - (nbands x 1) array of eigenenergies
%   eta - Broadening parameter for spectral functions
%   de  - Energy grid spacing for numerical integration
%   T   - (optional) Temperature in K, if not specified, take it from eta
% Output:
%   Imn - (nbands x nbands) is the -dfde * An * Am (spectral functions)

kB = 8.617333262145e-5;  % eV/K
if nargin == 3
    T = eta / kB; % use the broadening factor to determine temperature
end
% nbands = length(Enk);
kT = kB * T;
de = min(eta/10, kT/20);  % Finer of eta/10 or kT/20
de = max(de, 1e-6);       % Prevent excessively small spacing

% T=0 case: -df/dE = delta(E-EF)
if abs(kT) < 1e-14
    % Spectral functions at Fermi level
    denom_n = (EF - Enk).^2 + eta^2;
    A_n = (1/pi) * (eta ./ max(denom_n, 1e-16));
    Imn = A_n * A_n';
    return;
end

thermal_range = 6 * kT;     % Capture 99.9% of -df/dE
spectral_range = 5 * eta;   % Capture most of spectral function width

% Important energy regions: around EF and around each band energy
critical_points = [EF; Enk];
energy_margin = max(thermal_range, spectral_range);
energy_min = min(critical_points) - energy_margin;
energy_max = max(critical_points) + energy_margin;

% Ensure reasonable grid size (avoid too many points)
max_points = 10000;
current_points = (energy_max - energy_min) / de;
if current_points > max_points
    de = (energy_max - energy_min) / max_points;
end

energy_grid = (energy_min:de:energy_max)';
% npoints = length(energy_grid);

% VECTORIZED SPECTRAL FUNCTION CALCULATION
% Using implicit expansion for efficiency
epsilon_diff_sq = (energy_grid - Enk').^2;  % (nenergy_grids x nbands)
denom = epsilon_diff_sq + eta^2;
A = (1/pi) * (eta ./ denom)';              % (nbands x nenergy_grids)

% OPTIMIZED FERMI DERIVATIVE at EF
x = (energy_grid - EF) / kT;

% Single robust calculation avoiding intermediate large values
% -df/dE = e^{(E-EF)/kT} / [kT * (1 + e^{(E-EF)/kT})^2]
% This form is numerically stable for all x
x_pos = max(x, -50);  % Prevent underflow
exp_x = exp(x_pos);
dfde = exp_x ./ (kT * (1 + exp_x).^2);

% throw numerical artifacts away
dfde = max(dfde, 0);
dfde(isnan(dfde) | isinf(dfde)) = 0;

% EFFICIENT INTEGRATION USING MATRIX OPERATIONS
% Scale by sqrt(weights) and use single matrix multiplication
sqrt_w = sqrt(dfde');
A_scaled = A .* sqrt_w;
Imn = (A_scaled * A_scaled') * de;

% Final numerical cleanup
Imn = real(Imn);
Imn = (Imn + Imn') / 2;
Imn = max(Imn, 0);  % Physical constraint

end
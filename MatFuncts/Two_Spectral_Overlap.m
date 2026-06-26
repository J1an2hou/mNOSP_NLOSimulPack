function Imn = Two_Spectral_Overlap(Enk, eta, EF, T)
% EVALUATE_CONDUCTIVITY_INTEGRAL 
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
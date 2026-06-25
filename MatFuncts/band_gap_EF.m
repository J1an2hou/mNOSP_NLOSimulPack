function [band_gap, EF] = band_gap_EF(Enk, Ne_tot, T_kelvin,tol)
% This code determines the band gap and EF according to Enk
    % E: Matrix of size [nkpts x nbands]
    % Ne_tot: Scalar
    % T_kelvin: Temperature in Kelvin
    % tol: tolerance factor
    if nargin == 3
        tol = 1e-8;
    end
    % Constants
    kB = 8.617333262e-5; % Boltzmann constant in eV/K
    beta = 1 / (kB * T_kelvin);
    
    % Flatten energy levels for easier summation
    all_energies = Enk(:);
    N_k = size(Enk, 1); % Number of k-points
    nbands = size(Enk, 2);
    
    % Define the electron count function: f(mu) = (Sum of FD) / N_k
    % We normalize by N_k because each k-point represents a state
    calc_ne = @(mu) sum(1 ./ (exp((all_energies - mu) * beta) + 1)) / N_k;
    
    % 1. Determine Fermi Level (EF) using Bisection
    % Range: slightly outside the min/max of your spectrum
    low = min(all_energies) - 0.5;
    high = max(all_energies) + 0.5;
    
    for iter = 1:500
        mid = (low + high) / 2;
        if calc_ne(mid) < Ne_tot
            low = mid;
        else
            high = mid;
        end

        if (high - low) < tol
            break
        elseif iter == 500
            warning('Total iteration reached for EF evaluation')
        end
    end
EF = (low + high) / 2;
    
% Find band gap
% For each k-point, find the highest occupied and lowest unoccupied bands
max_occ_band = floor(Ne_tot);  % Assuming integer occupation per band
if max_occ_band >= nbands
    error('Total occupation exceeds available bands');
end

% Initialize arrays for valence band maximum (VBM) and conduction band minimum (CBM)
VBM = -inf;
CBM = inf;

% Loop over all k-points
for ik = 1:N_k
    % Sort energies at this k-point
    [~, idx] = sort(Enk(ik, :));
    
    % Find VBM (maximum energy among occupied bands)
    if max_occ_band > 0
        vbm_candidate = max(Enk(ik, idx(1:max_occ_band)));
        VBM = max(VBM, vbm_candidate);
    end
    
    % Find CBM (minimum energy among unoccupied bands)
    if max_occ_band < nbands
        cbm_candidate = min(Enk(ik, idx(max_occ_band+1:end)));
        CBM = min(CBM, cbm_candidate);
    end
end

% Calculate band gap
if isinf(VBM) || isinf(CBM)
    band_gap = 0;

else
    band_gap = CBM - VBM;
    
    % Check if system is metallic (considering temperature broadening)
    % If gap is very small compared to temperature, treat as metallic
    if band_gap < 10*(kB * T_kelvin)

        band_gap = 0;
    else

    end
end

end

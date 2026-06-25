function [Eg, gap_type, VBM, CBM] = band_gap(E, EF)
% CALCULATE_BAND_GAP Calculate band gap from eigenenergies and Fermi energy
%
% Input:
%   E - matrix of eigenenergies (nkpoints x nbands)
%   EF - Fermi energy
%
% Output:
%   band_gap - band gap value
%   gap_type - 'direct', 'indirect', or 'metal'
%   VBM - Valence Band Maximum energy
%   CBM - Conduction Band Minimum energy

    % Find valence and conduction bands
    % Valence bands: energies below or equal to EF
    % Conduction bands: energies above EF
    
    valence_bands = E(E <= EF);
    conduction_bands = E(E > EF);
    
    % Check if material is metallic
    if isempty(valence_bands) || isempty(conduction_bands)
        Eg = 0;
        gap_type = 'metal';
        VBM = EF;
        CBM = EF;
        return;
    end
    
    % Find VBM and CBM
    VBM = max(valence_bands);
    CBM = min(conduction_bands);
    
    % Calculate indirect gap (minimum overall gap)
    indirect_gap = CBM - VBM;
    
    % Calculate direct gaps at each k-point
    nkpoints = size(E, 1);
    direct_gaps = zeros(nkpoints, 1);
    
    for ik = 1:nkpoints
        kpoint_energies = E(ik, :);
        vbm_k = max(kpoint_energies(kpoint_energies <= EF));
        cbm_k = min(kpoint_energies(kpoint_energies > EF));
        
        if ~isempty(vbm_k) && ~isempty(cbm_k)
            direct_gaps(ik) = cbm_k - vbm_k;
        else
            direct_gaps(ik) = inf;
        end
    end
    
    direct_gap = min(direct_gaps);
    
    % Determine band gap type and value
    if indirect_gap <= 0
        Eg = 0;
        gap_type = 'metal';
    elseif abs(direct_gap - indirect_gap) < 1e-10  % Numerical tolerance
        Eg = direct_gap;
        gap_type = 'direct';
    else
        Eg = indirect_gap;
        gap_type = 'indirect';
    end
    
    % Display results
    % fprintf('Band gap: %.4f eV\n', Eg);
    % fprintf('Gap type: %s\n', gap_type);
    % fprintf('VBM: %.4f eV\n', VBM);
    % fprintf('CBM: %.4f eV\n', CBM);
end
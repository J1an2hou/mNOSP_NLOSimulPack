function Imn = Two_Spectral_Overlap_quadgk(Enk, eta, EF, T)
% Inm_precise  Hybrid precise evaluation of
%   I_nm = ∫ A(ε;Enk(n)) A(ε;Enk(m)) (-df/dε) dε
% where A is a normalized Lorentzian with width Gamma.
%
% Behavior:
%  - If T == 0: returns exact T->0 result I_nm = A(mu;Enn) * A(mu;Emm)
%  - For T > 0: uses analytic J_nm with df sampled at midpoint for
%    |Delta| > thr, and uses quadgk numeric integration for |Delta| <= thr.
%
% Inputs:
%   Enk  : [N x 1] band energies (eV)
%   Gamma: scalar broadening (eV)
%   mu   : chemical potential (eV)
%   T    : temperature (K)  (T==0 handled specially)
%
% Output:
%   I : [N x N] matrix of I_nm
%
% Notes:
%   - The analytic J_nm used is J = (2/pi)*Gamma / (Delta^2 + (2*Gamma)^2)
%   - The Fermi derivative used is -df/dE = (1/(4 kB T)) sech^2((E-mu)/(2 kB T))
%   - Numeric integration uses quadgk with tight tolerances.
%
% ensure column
Enk = Enk(:);
% N = length(Enk);

kB = 8.617333262145e-5;  % eV/K
if nargin == 3
    T = eta / kB; % use the broadening factor to determine temperature
end
% nbands = length(Enk);
kT = kB * T;

% T=0 case: -∂f/∂ε = δ(ε-EF)
if abs(kT) < 1e-14
    % Spectral functions at Fermi level
    denom_n = (EF - Enk).^2 + eta^2;
    A_n = (1/pi) * (eta ./ max(denom_n, 1e-16));
    Imn = A_n * A_n';
    return;
end

Delta = Enk - Enk.';            % N x N
E0    = 0.5*(Enk + Enk.');      % N x N midpoint

% Fermi-derivative sampled at midpoint (positive quantity)

X = (E0 - EF) ./ (2*kT);
dfE0 = (1./(4*kT)) .* sech(X).^2;   % N x N

% analytic spectral overlap J_nm (exact integral of A_n * A_m)
J = (2/pi) * ( eta ./ ( Delta.^2 + (2*eta).^2 ) );  % N x N

I_analytic = dfE0 .* J;   % N x N matrix (fast)

% --- choose threshold for numerical treatment (automatic) --------------
% threshold uses both Gamma and thermal width
thr = max(5*eta, 5*kT);   % you can adjust multiplier if desired

% find problematic pairs where |Delta| <= thr
mask_num = abs(Delta) <= thr;          % logical N x N
% we'll only integrate for upper triangle or diagonal and mirror
[i_idx, j_idx] = find(triu(mask_num));  % indices (i <= j)

% start with analytic result
Imn = I_analytic;

% If no pairs to integrate, return fast
if isempty(i_idx)
    % safety: ensure symmetry and real
    Imn = 0.5*(Imn + Imn.');
    Imn = real(Imn);
    return;
end

% --- prepare numeric integrator helpers ----------------------------------
% vectorized lorentzian and f-derivative
A_of_e = @(e, eps0) (1/pi) .* ( eta ./ ( (e - eps0).^2 + eta.^2 ) );
minus_df = @(e) (1./(4*kT)) .* sech( (e - EF)./(2*kT) ).^2;  % returns positive -df

% integration parameters (tight but reasonable)
relTol = 1e-8;
absTol = 1e-12;
maxWidthFactor = 20;   % integration half-width multiplier relative to max(Gamma,kBT)
% precompute a safe global window to avoid repeated large ranges:
globalHalfWidth = maxWidthFactor * max([eta, kT, thr]);
global_emin = EF - globalHalfWidth;
global_emax = EF + globalHalfWidth;

% --- perform numeric integrals for needed pairs (upper triangular only) ---
numPairs = length(i_idx);
for p = 1:numPairs
    n = i_idx(p);
    m = j_idx(p);
    En = Enk(n);
    Em = Enk(m);
    % choose local window centered at midpoint E0(n,m)
    center = 0.5*(En + Em);
    % local half-width should cover both spectral tails and thermal tails
    localHalf = max( max(8*eta, 8*kT), abs(center - EF) + 6*eta );
    emin = max(global_emin, center - localHalf);
    emax = min(global_emax, center + localHalf);
    % integrand function handle
    integrand = @(e) A_of_e(e, En) .* A_of_e(e, Em) .* minus_df(e);
    % integrate with quadgk
    try
        val = quadgk(integrand, emin, emax, 'RelTol', relTol, 'AbsTol', absTol);
    catch ME
        % fallback: expand window and retry once
        emin2 = EF - 50*max([eta,kT]);
        emax2 = EF + 50*max([eta,kT]);
        val = quadgk(integrand, emin2, emax2, 'RelTol', 1e-7, 'AbsTol',1e-11);
    end
    % assign symmetrically
    Imn(n,m) = val;
    Imn(m,n) = val;
end

% ensure symmetry and real small numerical imag parts removed
Imn = 0.5*(Imn + Imn.');
Imn = real(Imn);

end
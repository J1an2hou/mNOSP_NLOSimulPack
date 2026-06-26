function Pdir = Polarization_Wilson_layer(eigvecs, lproj, Nocc, Nk1, Nk2, Nk3, dir)
% ------------------------------------------------------------
% Berry-phase polarization along one reciprocal direction
%
% INPUT:
% eigvecs : (Nk x Norb x Norb)
% Nocc    : number of occupied bands
% Nk1,Nk2,Nk3 : k-mesh sizes (Nk1 fastest)
% dir     : 1, 2, or 3   (b1, b2, b3)
%
% OUTPUT:
% Pdir : polarization along dir (in units of e, modulo e)
% ------------------------------------------------------------

% ------------------------------------------------------------
% ------------------------------------------------------------
%
%                      !!! BIG WARNING !!!
% The layer projection may not be correct in the current scheme
% Until we fix it, don't trust it at all
%          !!! This is very very very important !!!
%
% ------------------------------------------------------------
% ------------------------------------------------------------

% You should manually check the bandgap to see if it is too small. In such
% a case the Wilson loop method is not applicable.

% In order to convert it to conventional polarization, you should multiply
% the lattice vector along the direction dir, and divide it by total volume
nwl = size(lproj, 3);

    idx = @(i1,i2,i3) (i3-1) + (i2-1)*Nk3 + (i1-1)*Nk2*Nk3 + 1;
    switch dir
        case 1
            Nkdir = Nk1;  Nperp = Nk2 * Nk3;
        case 2
            Nkdir = Nk2;  Nperp = Nk1 * Nk3;
        case 3
            Nkdir = Nk3;  Nperp = Nk1 * Nk2;
        otherwise
            error('dir must be 1, 2, or 3');
    end
    phase_sum = zeros(1,nwl);
    for iwl = 1:nwl
    for p = 1:Nperp
        W = eye(Nocc,Nocc);
        for l = 1:Nkdir
            [i1,i2,i3] = transverse_indices(p, l, dir, Nk2, Nk3);
            k  = idx(i1,i2,i3);
            switch dir
                case 1
                    kp = idx(mod(i1,Nk1)+1, i2, i3);
                case 2
                    kp = idx(i1, mod(i2,Nk2)+1, i3);
                case 3
                    kp = idx(i1, i2, mod(i3,Nk3)+1);
            end
        
            U  = squeeze(eigvecs(k,:,1:Nocc));
            Up = squeeze(eigvecs(kp,:,1:Nocc));
            M = U' * lproj(:,:,iwl) * Up; % if layer projection wanted, M = U' * lproj * Up;
            W = M * W;
        end
        phase_sum(iwl) = phase_sum(iwl) + imag(log(det(W)));
    end
    end
    Pdir = -phase_sum / (2*pi*Nperp);
    Pdir = mod(Pdir,1);  % polarization quantum, in unit of e (e = 1.6e-19 C)
    % This Pdir is mod 1. Hence, a 0.999 is the same as -0.001.

end


function [i1,i2,i3] = transverse_indices(p, l, dir, Nk2, Nk3)
% ------------------------------------------------------------
% k-loop index mapping
% Nk3 fastest, Nk1 slowest
% ------------------------------------------------------------

    switch dir
        case 1
            % Loop along i1
            i1 = l;
            i2 = mod(p-1, Nk2) + 1;
            i3 = floor((p-1)/Nk2) + 1;

        case 2
            % Loop along i2
            i2 = l;
            i3 = mod(p-1, Nk3) + 1;
            i1 = floor((p-1)/Nk3) + 1;

        case 3
            % Loop along i3
            i3 = l;
            i2 = mod(p-1, Nk2) + 1;
            i1 = floor((p-1)/Nk2) + 1;
    end
end


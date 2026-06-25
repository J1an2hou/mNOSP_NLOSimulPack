function [bopv_on_k,bspv_on_k] = ...
    bpv_vgauge_voltage(o,v,z,f,orbital_seq,nw,wavefunct,omega,ksi,direct)
ndir = size(direct,1);
nhv = length(omega);
bopv_on_k = zeros(ndir, nhv);
bspv_on_k = zeros(ndir, nhv);
%jorb = zeros(ndir,nw);
%jspin = zeros(ndir,nw);
 
for idir = 1:ndir
a = direct(idir,1); % here a should be along non-periodic direction
b = direct(idir,2);
c = direct(idir,3);
d = direct(idir,4);
if a ~= 3
    warning('This code is for non-periodic voltage.')
end
va = -z; % the negative sign comes from <dipole> = -e*z
%ra = r(:,:,a);
od = calc_orb(nw,wavefunct,orbital_seq,d);
sd = calc_sp(nw,wavefunct,d);
% ====== Here two different definitions can be used ======

jod = (od*va+va*od)/2;                  % This is conventional definition
%jod = 1i*(ra*(omn.*od)+(omn.*od)*ra)/2; % This is the added torque term
%jod = 1i*omn .* (ra*od+od*ra)/2;  % This is d(r*L)/dt, anticommutation

jsd = (sd*va+va*sd)/2;                  % This is conventional definition
%jsd = 1i*(ra*(omn.*sd)+(omn.*sd)*ra)/2; % This is the added torque term
%jsd = 1i*omn .* (ra*sd+sd*ra)/2;  % This is d(r*S)/dt, anticommutation

% ========================================================
%jorb(idir,:) = diag(real(jod));
%jspin(idir,:) = diag(real(jsd));
facC = v(:,:,c);

facBorb = jod ./ (o+1i*ksi);
facBspin = jsd ./ (o+1i*ksi);

BCCBorb = facBorb*facC-facC*facBorb;
BCCBspin = facBspin*facC-facC*facBspin;
for iomega = 1:length(omega)
    o_in = omega(iomega);
    facA = f .* v(:,:,b) ./(-o-o_in+1i*ksi);
    bopv_on_k(idir, iomega) = trace(facA * BCCBorb);
    bspv_on_k(idir, iomega) = trace(facA * BCCBspin);
end
    if b == c   % linearly polarized light
        bopv_on_k(idir,:) = real(bopv_on_k(idir,:));
        bspv_on_k(idir,:) = real(bspv_on_k(idir,:));
    else        % circularly polarized light
        bopv_on_k(idir,:) = imag(bopv_on_k(idir,:));
        bspv_on_k(idir,:) = imag(bspv_on_k(idir,:));
    end
end
end
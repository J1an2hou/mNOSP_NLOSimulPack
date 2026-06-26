function [intersp_on_k,interorb_on_k,intrasp_on_k,intraorb_on_k] ...
    = LinearEdel(vmn,omn,lpmn,fmn,wavefunct,orbital_seq,nw,hv,direct,dfde,ksi)
nhv = length(hv);
ndir = size(direct,1);
intersp_on_k = zeros(ndir, nhv);
interorb_on_k = zeros(ndir, nhv);
intrasp_on_k = zeros(ndir, nhv);
intraorb_on_k = zeros(ndir, nhv);

for idir = 1:ndir
    a = direct(idir,1);
    b = direct(idir,2);
    vnma = vmn(:,:,a).';
    orb_i = calc_orb(nw,wavefunct,orbital_seq,b)*lpmn;
    sp_i = calc_sp(nw,wavefunct,b)*lpmn;
    fac1spinter = fmn .* omn ./ (omn.^2+ksi^2) .* sp_i;
    fac1orbinter = fmn .* omn ./ (omn.^2+ksi^2) .* orb_i;
    for ihv = 1:nhv
        o = hv(ihv);
        fac2inter = vnma ./ (o - omn' + 1i*ksi);
        intersp_on_k(idir,ihv) = trace(fac1spinter*fac2inter);
        interorb_on_k(idir,ihv) = trace(fac1orbinter*fac2inter);
        intrasp_on_k(idir,ihv) = sum(dfde.*diag(sp_i).*diag(vnma)/(o+1i*ksi));
        intraorb_on_k(idir,ihv) = sum(dfde.*diag(orb_i).*diag(vnma)/(o+1i*ksi));
    end
end

end
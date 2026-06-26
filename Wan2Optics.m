% This code is read and analyze wannier90_hr.dat file
% Everything in eV or Angstrom

clear;

%% Initial
nkx = 301;
nky = 301;
nkz = 1;
spin_d = 1;   % spin degeneracy: 1 - spin polarized or soc; 2 - unpolarized
ksi = 0.02;  % damping term, in eV
eps0 = 1/181;    % vacuum permittivity, eps_0=e^2/(2\alpha*hc)
hv = linspace(0.1, 6, 1500); % optical frequency;
nocc = 20; % occupied band numbers
EF = -0.6902;

thickness = 6; % in Angstrom

%% Read Files
[R,HR,wtR,wc,nw,nR] = getwanhr('wannier90');

[alatt, blatt] = readlattice('POSCAR');
Vlatt = det(alatt);

%% Main
nhv = length(hv);
kx = linspace(-0.5, 0.5, nkx+1);
ky = linspace(-0.5, 0.5, nky+1);
kz = linspace(-0.5, 0.5, nkz+1);
dkx = 1/nkx; dky = 1/nky; dkz = 1/nkz;
kx = kx+dkx/2; ky = ky+dky/2; kz = kz+dkz/2;
[kx, ky, kz] = meshgrid(kx(1:nkx), ky(1:nky), kz(1:nkz));
kx = permute(kx,[2,1,3]);
ky = permute(ky,[2,1,3]);
kz = permute(kz,[2,1,3]);
bkx = zeros(nkx,nky,nkz);
bky = zeros(nkx,nky,nkz);
bkz = zeros(nkx,nky,nkz);
nkpts = nkx*nky*nkz;
kweight = 1/nkpts;
kptcell = cell(nkpts,1);
for ikpt = 1:nkpts
    ikx = ceil(ikpt/(nky*nkz)); % ikz is the fastest index
    iky = ceil((ikpt-(ikx-1)*nky*nkz)/nkz);
    ikz = ceil(ikpt-(ikx-1)*nky*nkz-(iky-1)*nkz);
    kptcell{ikpt} = [kx(ikx,iky,ikz), ky(ikx,iky,ikz), kz(ikx,iky,ikz)];
end
f = zeros(nw, 1);

chixx = zeros(1, nhv);
chiyy = zeros(1, nhv);
chixx_on_k = zeros(nkpts, nhv);
chiyy_on_k = zeros(nkpts, nhv);
chixx_FS = zeros(1, nhv);
chiyy_FS = zeros(1, nhv);
chixx_FS_on_k = zeros(nkpts, nhv);
chiyy_FS_on_k = zeros(nkpts, nhv);
jdos_on_k = zeros(nkpts, nhv);
jdos = zeros(1,length(hv));

parpool(4)
parfor ikpt = 1:nkpts
    disp(['ikpt = ',num2str(ikpt),' in total ',num2str(nkpts),' points']);
    kpt = kptcell{ikpt};
    % build Hamiltonian
    [Hk,dHdkx,dHdky,dHdkz] = buildHk(R,HR,wtR,wc,nw,nR,kpt,alatt,blatt);
      
    % diagonalize
    [V, D] = eig(Hk);
    [d,ind] = sort(real(diag(D)));
    energy = d;
    wavefunct = V(:,ind);
    vmn = calc_vmn(wavefunct, dHdkx, dHdky, dHdkz, nw);
    [rmn, omn] = calc_rmn(vmn,energy,nw,ksi);
    f = 1./(1+exp((energy-EF)/ksi));
    fmn = f*ones(1,nw)-ones(nw,1)*f';
    % calcuate the contribution to bpv from this k point
    [cxx,cyy,czz] = chi_interband(omn,vmn,fmn,hv,ksi);
    [cxx_FS,cyy_FS,czz_FS] = chi_intraband(vmn,fmn,hv,ksi);
    % calculate joint DOS
    % jDOS
    jd = zeros(1, length(hv));
    for ihv = 1:length(hv)
        domn = df(omn-hv(ihv),ksi);
        jd(ihv) = trace(fmn*domn);
    end
    jdos = jdos + kweight *jd /Vlatt*spin_d;
    chixx_on_k(ikpt, :) = cxx + cxx_FS; % interband plus intraband
    chiyy_on_k(ikpt, :) = cyy + cyy_FS;
    jdos_on_k(ikpt, :) = jd;
    chixx = chixx + kweight * (cxx+cxx_FS) /Vlatt/eps0*spin_d; % interband plus intraband
    chiyy = chiyy + kweight * (cyy+cyy_FS) /Vlatt/eps0*spin_d;
end

chixr = real(chixx) * alatt(3,3) / thickness;
chixi = 1-exp(-hv.*imag(chixx) * alatt(3,3)/12400);

chiyr = real(chiyy) * alatt(3,3) / thickness;
chiyi = 1-exp(-hv.*imag(chiyy) * alatt(3,3)/12400);

jdos_on_k = jdos_on_k * alatt(3,3) / thickness;
jdos = jdos * alatt(3,3) / thickness;

semilogx(hv,chixi*100,'b','LineWidth',1.5)
hold on
semilogx(hv,chiyi*100,'g','LineWidth',1.5)
xlabel('Photon energy ($\mathrm{eV}$)','interpreter','latex', ...
    'fontname','times new roman','fontsize',20)
ylabel('Absorption (\%)', ...
    'fontname','times new roman','fontsize',20)
legend('$x$-LPL','$y$-LPL','interpreter','latex', ...
    'fontname','times new roman','fontsize',20)
ax = gca;
ax.FontName='Times New Roman';
ax.LineWidth=1.5;
ax.FontSize=20;
%save data-101201 chixx_on_k chiyy_on_k jdos_on_k ...
%    chixr chixi chiyr chiyi jdos omega_all bkx bky
%hold off
%[X,Y]=meshgrid(kx,ky);
%colormap('parula');
%Zx=zeros(length(kx),length(ky));
%Zz=zeros(length(kx),length(ky));
%eg=zeros(length(kx),length(ky));
%for i=1:length(kx)
%for j=1:length(ky)
%Zx(i,j)=real(chix(1,i,j,1));
%Zz(i,j)=real(chiz(1,i,j,1));
%eg(i,j)=real(energy(10,i,j,1)-energy(9,i,j,1));
%end
%end
%s=surf(X,Y,real(eg),real(Zx),'EdgeColor','interp',...
%  'FaceColor','interp','FaceAlpha',0.75,'FaceLighting','gouraud');
%s.LineStyle='none';
%h = light;
%h.Style='infinite';
%h.Position=[-1 -1 -1];
%xlabel('\it k_{x}','FontSize',18,'FontName','Times New Roman')
%ylabel('\it k_{y}','FontSize',18,'FontName','Times New Roman')
%xticks(-0.4:0.2:0.4)
%yticks(-0.4:0.2:0.4)
%zlabel('Bandgap (eV)','FontSize',18,'FontName','Times New Roman')
%ax = gca;
%ax.FontName='Times New Roman';
%set(gca,'Fontname','Times New Roman');
%set(gca,'FontSize',18);
%c = colorbar;
%c.Label.String = '\Re\chi_{\it xx}(\it k\rm)';


%% Functions
function [chixx_k,chiyy_k,chizz_k] = chi_interband(omn,vmn,fmn,omega_all,omega0)
chixx_k = zeros(1, length(omega_all));
chiyy_k = zeros(1, length(omega_all));
chizz_k = zeros(1, length(omega_all));
Amn1 = -1i*vmn(:,:,1).*omn./(omn.^2+omega0^2);
Amn2 = -1i*vmn(:,:,2).*omn./(omn.^2+omega0^2);
Amn3 = -1i*vmn(:,:,3).*omn./(omn.^2+omega0^2);
for ihv = 1:length(omega_all)
   o = omega_all(ihv);
   fac = omn-o-1i*omega0;
   chixx_k(ihv) = trace((Amn1./fac)*(fmn.*Amn1));
   chiyy_k(ihv) = trace((Amn2./fac)*(fmn.*Amn2));
   chizz_k(ihv) = trace((Amn3./fac)*(fmn.*Amn3));
end

end

function [chixx_intra, chiyy_intra, chizz_intra] = chi_intraband(vmn,energy,omega_all,omega0,EF)
% chi = e^2/eps0 * int_k (v_nn^a * v_nn^a * delta(En-EF)/omega/(omega+i*eta) )
chixx_intra = zeros(1, length(omega_all));
chiyy_intra = zeros(1, length(omega_all));
chizz_intra = zeros(1, length(omega_all));
vx = diag(real(vmn(:,:,1)));
vy = diag(real(vmn(:,:,2)));
vz = diag(real(vmn(:,:,3)));
for ihv = 1:length(omega_all)
    o = omega_all(ihv);
    facx = sum((vx.^2) .* df(energy-EF,omega0));
    facy = sum((vy.^2) .* df(energy-EF,omega0));
    facz = sum((vz.^2) .* df(energy-EF,omega0));
    fac = (o+1i*omega0)*o;
    chixx_intra(ihv) = facx * fac / (fac^2+omega0^2);
    chiyy_intra(ihv) = facy * fac / (fac^2+omega0^2);
    chizz_intra(ihv) = facz * fac / (fac^2+omega0^2);
end


end

function vmn = calc_vmn(wavefunct, dHdkx, dHdky, dHdkz, nw)
vmn = zeros(nw, nw, 3);

vmn(:,:,1) = wavefunct' * dHdkx * wavefunct;
vmn(:,:,2) = wavefunct' * dHdky * wavefunct;
vmn(:,:,3) = wavefunct' * dHdkz * wavefunct;

end

function [rmn, omn] = calc_rmn(vmn,energy,nw,eta)
rmn = zeros(nw,nw,3);
omn = energy*ones(1,nw) - ones(nw,1)*energy';
for i = 1:nw
    for j = 1:nw
        
        if i ~= j
            rmn(i,j,:) = -1i*vmn(i,j,:)*omn(i,j)/(omn(i,j)^2+eta^2);
        end
    end
end
end

function deltafunct = df(x,ksi)
deltafunct = ksi./(x.^2+ksi^2)/pi; % Lorentz type
% deltafunct = 1/ksi/sqrt(pi)*exp(-x.^2/ksi^2); % Gaussian type
end

function [alatt, blatt] = readlattice(filename)
file = fopen(filename);
fgetl(file); 
scale = sscanf( fgetl(file), '%f' );
alatt = zeros(3,3);
for i = 1:3
  alatt(i,:) = sscanf( fgetl(file), '%f' );
end
alatt = scale * alatt;
V = det(alatt);

blatt = zeros(3,3);
blatt(1,:) = 2*pi*cross(alatt(2,:),alatt(3,:))/V;
blatt(2,:) = 2*pi*cross(alatt(3,:),alatt(1,:))/V;
blatt(3,:) = 2*pi*cross(alatt(1,:),alatt(2,:))/V;

end


function [R,HR,wtR,wc,nw,nR] = getwanhr(seedname)
file = fopen(strcat([seedname, '_hr.dat']), 'r');
fgetl(file);
line = fgetl(file);
nw = sscanf(line, '%d');
line = fgetl(file);
nR = sscanf(line, '%d');

Rdeg = [];
Rlines = ceil(nR/15); % per the format of wannier90_hr, 15 entries per line
for iline = 1:Rlines
	line = fgetl(file);
  Rdeg = [Rdeg, sscanf(line, '%d')'];
end
wtR = 1./Rdeg;

C = textscan(file,'%d %d %d %d %d %n %n');
R1 = C{1};
R2 = C{2};
R3 = C{3};
Hr1 = C{6};
Hr2 = C{7};
% nlines = length(R1);
% nw = sqrt(nlines/nR);
fclose(file);
R = zeros(nR,3);
HR = zeros(nR,nw,nw);
for i = 1:nR
    R(i,1) = R1((i-1)*nw^2+1);
    R(i,2) = R2((i-1)*nw^2+1);
    R(i,3) = R3((i-1)*nw^2+1);
    for j = 1:nw
        for k = 1:nw
            HR(i,j,k) = Hr1((i-1)*nw^2+(k-1)*nw+j) ...
                +1i*Hr2((i-1)*nw^2+(k-1)*nw+j);
        end
    end
    
end
% symmetrize HR
for i = 1:nR
  HR(i,:,:) = 0.5*( HR(i,:,:) + conj(permute(HR(nR+1-i,:,:), [1,3,2])) );
end
% wc = zeros(nw,3);

fileID = fopen(strcat([seedname,'_centres.xyz']),'r');
fgetl(fileID);
fgetl(fileID);
C = textscan(fileID,'%*s %f %f %f');
wc1 = C{1};
wc2 = C{2};
wc3 = C{3};
fclose(fileID);
wc = zeros(nw,3);
for i = 1:nw
    wc(i,1) = wc1(i);
    wc(i,2) = wc2(i);
    wc(i,3) = wc3(i);
end

end

function [Hk,dHdkx,dHdky,dHdkz] = buildHk(R,HR,wtR,wc,nw,nR,kpt,alatt,blatt)
% Note that here kpoint coordinate is in direct mode (in unit of reciprocal
% space
Hk = zeros(nw,nw);
dHdkx = zeros(nw,nw);
dHdky = zeros(nw,nw);
dHdkz = zeros(nw,nw);
bk = kpt*blatt;
for iR = 1:nR
    wcdiff1 = wc(:,1)*ones(1,nw)-ones(nw,1)*wc(:,1)';
    wcdiff2 = wc(:,2)*ones(1,nw)-ones(nw,1)*wc(:,2)';
    wcdiff3 = wc(:,3)*ones(1,nw)-ones(nw,1)*wc(:,3)';
    aR = R(iR,:)*alatt;
    aR1 = aR(1) - wcdiff1;
    aR2 = aR(2) - wcdiff2;
    aR3 = aR(3) - wcdiff3;
    kdotR = aR1 * bk(1) + aR2 * bk(2) + aR3 * bk(3);
    hr = reshape(HR(iR,:,:),nw,nw);
    Hk = Hk + wtR(iR) * hr .* exp(1i*kdotR);
    dHdkx = dHdkx + 1i*wtR(iR)*aR1 .* hr .* exp(1i*kdotR);
    dHdky = dHdky + 1i*wtR(iR)*aR2 .* hr .* exp(1i*kdotR);
    dHdkz = dHdkz + 1i*wtR(iR)*aR3 .* hr .* exp(1i*kdotR);
end
 
rdn = -6;
Hk = roundn(Hk,rdn);
dHdkx = roundn(dHdkx,rdn);
dHdky = roundn(dHdky,rdn);
dHdkz = roundn(dHdkz,rdn);
end
 



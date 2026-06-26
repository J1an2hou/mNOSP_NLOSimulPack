% This code is read and analyze wannier90_hr.dat file
% Everything in eV or Angstrom

clear all;

% Initial parameters
nocc = 2; % occupied band numbers
EF = -3:0.01:-2;
dkx = 0.01;
dky = 0.01;
dkz = 0.45;
spin_d = 1;   % spin degeneracy: 1 - spin polarized or soc; 2 - unpolarized
ksi = 0.01;  % dampling term, in eV
eps0 = 1.0/181;    % vacuum permittivity, eps_0=e^2/(2\alpha*hc)
alatt = [4 0.0 0.0; ... % real lattice
         -2 3.4641 0.0; ...
         0.0 0.0 15.0];
kx = -0.5+dkx/2:dkx:0.5-dkx/2;
ky = -0.5+dky/2:dky:0.5-dky/2;
kz = -0.5+dkz/2:dkz:0.5-dkz/2;


V = dot(alatt(1,:),cross(alatt(2,:),alatt(3,:)));
blatt = zeros(3,3);
blatt(1,:) = 2*pi*cross(alatt(2,:),alatt(3,:))/V;
blatt(2,:) = 2*pi*cross(alatt(3,:),alatt(1,:))/V;
blatt(3,:) = 2*pi*cross(alatt(1,:),alatt(2,:))/V;
b_s = blatt(1,1)*blatt(2,2) - blatt(1,2)*blatt(2,1);

[R,HR,wtR,wc,nw,nR] = getwanhr('wannier90-');
ikpt = 0;
Hk = zeros(nw,nw,length(kx),length(ky),length(kz));
dHdkx = zeros(nw,nw,length(kx),length(ky),length(kz));
dHdky = zeros(nw,nw,length(kx),length(ky),length(kz));
dHdkz = zeros(nw,nw,length(kx),length(ky),length(kz));
energy = zeros(nw,length(kx),length(ky),length(kz));
wavefunct = zeros(nw,nw,length(kx),length(ky),length(kz));
for ikx = 1:length(kx)
    for iky = 1:length(ky)
        for ikz = 1:length(kz)
            ikpt = ikpt+1;
            kpt = [kx(ikx),ky(iky),kz(ikz)];
            bkpt = kpt*blatt;
            [Hk(:,:,ikx,iky,ikz),dHdkx(:,:,ikx,iky,ikz), ...
                dHdky(:,:,ikx,iky,ikz),dHdkz(:,:,ikx,iky,ikz)] ...
                = buildHk(R,HR,wtR,wc,nw,nR,kpt,alatt);
            [Vect,D] = eig(Hk(:,:,ikx,iky,ikz));
            [d,ind] = sort(diag(D));
            Ds = D(ind,ind);
            for i = 1:nw
                energy(i,ikx,iky,ikz) = Ds(i,i);
            end
            wavefunct(:,:,ikx,iky,ikz) = Vect(:,ind);
        end
    end
end


[Amnx,Amny,Amnz,Asmnx,Asmny,Asmnz] ...
    = udhdku(wavefunct,dHdkx,dHdky,dHdkz,energy,kx,ky,kz,nw);
%[Amnx,Amny,Amnz,Asmnx,Asmny,Asmnz] ...
%    = udu(wavefunct,kx,ky,kz,dkx,dky,dkz,blatt,nw);

[vx,vy,vz] = dEdk(energy,kx,ky,kz,blatt,nw);

%[X,Y]=meshgrid(kx,ky);
%xx=zeros(length(kx),length(ky));
%yy=xx;
%for i=1:4:length(kx)
%for j=1:4:length(ky)
%xx(i,j)=vx(2,i,j,1);
%yy(i,j)=vx(2,i,j,1);
%end
%end
%quiver(X,Y,xx,yy)

% Now one can integrate everything
hv = 0:0.02:5;
[chixx,chiyy,chizz,chixy,chiyz,chizx] ...
    = chi(Amnx,Amny,Amnz,energy,hv,kx,ky,kz,nw,nocc,ksi,V,spin_d);

gmh = zeros(length(EF),1);
sgmh = gmh;
for imu = 1:length(EF)
[BCx,BCy,BCz,sBCx,sBCy,sBCz] ...
    = Berrycurv(Amnx,Amny,Amnz,Asmnx,Asmny,Asmnz ...
    ,kx,ky,kz,nw,nocc,energy,EF(imu));

[G_MH,sG_MH] = MagnusHall(BCz,sBCz,vx,kx,ky,nw,energy,EF(imu));
gmh(imu) = G_MH;
sgmh(imu) = sG_MH;
end

plot(EF,gmh,'r')
hold on
plot(EF,sgmh,'b','LineWidth',2)

chiallxx = zeros(length(hv),1);
chiallyy = zeros(length(hv),1);
chiallzz = zeros(length(hv),1);

for ihv = 1:length(hv)
    chiallxx(ihv) = chiallxx(ihv)+sum(sum(sum(chixx(ihv,:,:,:))));
    chiallyy(ihv) = chiallyy(ihv)+sum(sum(sum(chiyy(ihv,:,:,:))));
    chiallzz(ihv) = chiallzz(ihv)+sum(sum(sum(chizz(ihv,:,:,:))));
end



%semilogx(hv,real(chiallxx),'b','LineWidth',2)
%hold on
%semilogx(hv,imag(chiallxx),'g','LineWidth',1)


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
%eg(i,j)=real(energy(nocc+1,i,j,1)-energy(nocc,i,j,1));
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

function [vx,vy,vz] = dEdk(energy,kx,ky,kz,blatt,nw)
vx = zeros(nw,length(kx),length(ky),length(kz));
vy = zeros(nw,length(kx),length(ky),length(kz));
vz = zeros(nw,length(kx),length(ky),length(kz));
E = zeros(length(kx),length(ky),length(kz));
dkx = kx(2)-kx(1);
dky = ky(2)-ky(1);
dkz = kz(2)-kz(1);
n1 = blatt(1,:);
n2 = blatt(2,:);
n3 = blatt(3,:);
n1 = n1/dot(n1,n1);
n2 = n2/dot(n2,n2);
n3 = n3/dot(n3,n3);
b = [n1;n2;n3];
for iband = 1:nw
    for ikx = 1:length(kx)
    for iky = 1:length(ky)
    for ikz = 1:length(kz)
    E(ikx,iky,ikz) = real(energy(iband,ikx,iky,ikz));
    end
    end
    end
    [v1,v2,v3] = gradient(E,dkx,dky,dkz);
    for ikx = 1:length(kx)
        for iky = 1:length(ky)
            for ikz = 1:length(kz)
                A = b \ [v1(ikx,iky,ikz); ...
                         v2(ikx,iky,ikz); ...
                         v3(ikx,iky,ikz)];
                vx(iband,ikx,iky,ikz) = A(1);
                vy(iband,ikx,iky,ikz) = A(2);
                vz(iband,ikx,iky,ikz) = A(3);
            end
        end
    end
end

end

function [G_MH,sG_MH] = MagnusHall(BCz,sBCz,vx,kx,ky,nw,energy,EF)
% The Magnus Hall effect is evaluated according to PRL 123, 216802 (2019)
% The result is in unit of e^2/h * \Delta U / 2pi
% Note that the system must be 2D, so kz = 0.
G_MH = 0;
sG_MH = 0;
ikz = 1;
for iband = 1:nw
for ikx = 1:length(kx)
    for iky = 1:length(ky)
        if vx(iband,ikx,iky,ikz) > 0
            G_MH = G_MH + BCz(ikx,iky,ikz) ...
                *deltafunct(real(energy(iband,ikx,iky,ikz))-EF) ...
                /(length(kx)*length(ky));
            sG_MH = sG_MH + sBCz(ikx,iky,ikz) ...
                *deltafunct(real(energy(iband,ikx,iky,ikz))-EF) ...
                /(length(kx)*length(ky));
        end
    end
end
end

end

function y = deltafunct(x)
eps = 0.01;
y = exp(-x^2/4/eps)/2/sqrt(pi*eps);
end

function [BCx,BCy,BCz,sBCx,sBCy,sBCz] ...
    = Berrycurv(Amnx,Amny,Amnz,Asmnx,Asmny,Asmnz ...
    ,kx,ky,kz,nw,nocc,energy,EF)

BCx = zeros(length(kx),length(ky),length(kz));
BCy = zeros(length(kx),length(ky),length(kz));
BCz = zeros(length(kx),length(ky),length(kz));
sBCx = zeros(length(kx),length(ky),length(kz));
sBCy = zeros(length(kx),length(ky),length(kz));
sBCz = zeros(length(kx),length(ky),length(kz));
for ikx = 1:length(kx)
for iky = 1:length(ky)
for ikz = 1:length(kz)
    for iband = 1:nw
        if energy(iband,ikx,iky,ikz) > EF
            continue
        end
    for jband = 1:nw
        if energy(jband,ikx,iky,ikz) < EF
            continue
        end
        BCx(ikx,iky,ikz) = BCx(ikx,iky,ikz) + ...
            2*imag(Amny(iband,jband,ikx,iky,ikz) * ...
            Amnz(jband,iband,ikx,iky,ikz));
        BCy(ikx,iky,ikz) = BCy(ikx,iky,ikz) + ...
            2*imag(Amnz(iband,jband,ikx,iky,ikz) * ...
            Amnx(jband,iband,ikx,iky,ikz));
        BCz(ikx,iky,ikz) = BCz(ikx,iky,ikz) + ...
            2*imag(Amnx(iband,jband,ikx,iky,ikz) * ...
            Amny(jband,iband,ikx,iky,ikz));
        sBCx(ikx,iky,ikz) = sBCx(ikx,iky,ikz) + ...
            2*imag(Asmny(iband,jband,ikx,iky,ikz) * ...
            Amny(jband,iband,ikx,iky,ikz));
        sBCy(ikx,iky,ikz) = sBCy(ikx,iky,ikz) + ...
            2*imag(Asmnz(iband,jband,ikx,iky,ikz) * ...
            Amnz(jband,iband,ikx,iky,ikz));
        sBCz(ikx,iky,ikz) = sBCz(ikx,iky,ikz) + ...
            2*imag(Asmnx(iband,jband,ikx,iky,ikz) * ...
            Amnx(jband,iband,ikx,iky,ikz));
    end
    end
end
end
end

end


function [chixx,chiyy,chizz,chixy,chiyz,chizx] ...
    = chi(Amnx,Amny,Amnz,energy,hv,kx,ky,kz,nw,nocc,ksi,V,spin_d)
eps0 = 1.0/181;
chixx = zeros(length(hv),length(kx),length(ky),length(kz));
chiyy = zeros(length(hv),length(kx),length(ky),length(kz));
chizz = zeros(length(hv),length(kx),length(ky),length(kz));
chixy = zeros(length(hv),length(kx),length(ky),length(kz));
chiyz = zeros(length(hv),length(kx),length(ky),length(kz));
chizx = zeros(length(hv),length(kx),length(ky),length(kz));
% Note chiyx = conj(chixy), chizy = conj(chiyz), and chixz = conj(chizx)
for ihv = 1:length(hv)
    for iband = 1:nocc
    for jband = nocc+1:nw
        for ikx = 1:length(kx)
        for iky = 1:length(ky)
        for ikz = 1:length(kz)
            chixx(ihv,ikx,iky,ikz) = chixx(ihv,ikx,iky,ikz) ...
               +(Amnx(iband,jband,ikx,iky,ikz) ...
               *Amnx(jband,iband,ikx,iky,ikz)) ...
               ./(energy(jband,ikx,iky,ikz)-energy(iband,ikx,iky,ikz) ...
               -hv(ihv)-1i*ksi) ...
               /(length(kx)*length(ky)*length(kz))/eps0/V*spin_d;
           chiyy(ihv,ikx,iky,ikz) = chiyy(ihv,ikx,iky,ikz) ...
               +(Amny(iband,jband,ikx,iky,ikz) ...
               *Amny(jband,iband,ikx,iky,ikz)) ...
               ./(energy(jband,ikx,iky,ikz)-energy(iband,ikx,iky,ikz) ...
               -hv(ihv)-1i*ksi) ...
               /(length(kx)*length(ky)*length(kz))/eps0/V*spin_d;
           chizz(ihv,ikx,iky,ikz) = chizz(ihv,ikx,iky,ikz) ...
               +(Amnz(iband,jband,ikx,iky,ikz) ...
               *Amnz(jband,iband,ikx,iky,ikz)) ...
               ./(energy(jband,ikx,iky,ikz)-energy(iband,ikx,iky,ikz) ...
               -hv(ihv)-1i*ksi) ...
               /(length(kx)*length(ky)*length(kz))/eps0/V*spin_d;
           chixy(ihv,ikx,iky,ikz) = chixy(ihv,ikx,iky,ikz) ...
               +(Amnx(iband,jband,ikx,iky,ikz) ...
               *Amny(jband,iband,ikx,iky,ikz)) ...
               ./(energy(jband,ikx,iky,ikz)-energy(iband,ikx,iky,ikz) ...
               -hv(ihv)-1i*ksi) ...
               /(length(kx)*length(ky)*length(kz))/eps0/V*spin_d;
           chiyz(ihv,ikx,iky,ikz) = chizz(ihv,ikx,iky,ikz) ...
               +(Amny(iband,jband,ikx,iky,ikz) ...
               *Amnz(jband,iband,ikx,iky,ikz)) ...
               ./(energy(jband,ikx,iky,ikz)-energy(iband,ikx,iky,ikz) ...
               -hv(ihv)-1i*ksi) ...
               /(length(kx)*length(ky)*length(kz))/eps0/V*spin_d;
           chizx(ihv,ikx,iky,ikz) = chizz(ihv,ikx,iky,ikz) ...
               +(Amnz(iband,jband,ikx,iky,ikz) ...
               *Amnx(jband,iband,ikx,iky,ikz)) ...
               ./(energy(jband,ikx,iky,ikz)-energy(iband,ikx,iky,ikz) ...
               -hv(ihv)-1i*ksi) ...
               /(length(kx)*length(ky)*length(kz))/eps0/V*spin_d;        
        end
        end
        end
    end
    end
end
end

function [Amnx,Amny,Amnz,Asmnx,Asmny,Asmnz] ...
    = udhdku(wavefunct,dHdkx,dHdky,dHdkz,energy,kx,ky,kz,nw)
% We use i * < u_n | dH/dk | u_m > / ( E_m - E_n ) to evaluate interband 
% Berry connection (E_n != E_m)
Amnx = zeros(nw,nw,length(kx),length(ky),length(kz));
Amny = zeros(nw,nw,length(kx),length(ky),length(kz));
Amnz = zeros(nw,nw,length(kx),length(ky),length(kz));
% We define Asmn_x = i * < u_n | {dHdk_y,s_z} | u_m > / (E_m - E_n) to
% evaluate spin Berry connection
Asmnx = zeros(nw,nw,length(kx),length(ky),length(kz));
Asmny = zeros(nw,nw,length(kx),length(ky),length(kz));
Asmnz = zeros(nw,nw,length(kx),length(ky),length(kz));
sz = zeros(nw,nw);
sx = zeros(nw,nw);
sy = zeros(nw,nw);
if mod(nw,2) == 0
    %for ib = 1:2:nw-1
    %    sz(ib,ib) = 1/2;
    %    sz(ib+1,ib+1) = -1/2;
    %    sx(ib,ib+1) = 1/2;
    %    sx(ib+1,ib) = 1/2;
    %    sy(ib,ib+1) = -1i/2;
    %    sy(ib+1,ib) = 1i/2;
    %end
    for ib = 1:nw/2
        sz(ib,ib) = 1/2;
        sz(ib+nw/2,ib+nw/2) = -1/2;
        sx(ib,ib+nw/2) = 1/2;
        sx(ib+nw/2,ib) = 1/2;
        sy(ib,ib+nw/2) = -1i/2;
        sy(ib+nw/2,ib) = 1i/2;
    end
end
en_deg = 1e-5; % energy degenerate tolerance
for ikx = 1:length(kx)
    for iky = 1:length(ky)
        for ikz = 1:length(kz)
            Amnx(:,:,ikx,iky,ikz) = ...
                wavefunct(:,:,ikx,iky,ikz)'* ...
                dHdkx(:,:,ikx,iky,ikz) ...
               *wavefunct(:,:,ikx,iky,ikz);
            Amny(:,:,ikx,iky,ikz) = ...
               wavefunct(:,:,ikx,iky,ikz)'* ...
               dHdky(:,:,ikx,iky,ikz) ...
               *wavefunct(:,:,ikx,iky,ikz);
            Amnz(:,:,ikx,iky,ikz) = ...
                wavefunct(:,:,ikx,iky,ikz)'* ...
                dHdkz(:,:,ikx,iky,ikz) ...
                *wavefunct(:,:,ikx,iky,ikz);
            Asmnx(:,:,ikx,iky,ikz) = ...
                wavefunct(:,:,ikx,iky,ikz)'* ...
                (dHdky(:,:,ikx,iky,ikz)*sz + ...
                sz*dHdky(:,:,ikx,iky,ikz)) * ...
                wavefunct(:,:,ikx,iky,ikz)/2;
            Asmny(:,:,ikx,iky,ikz) = ...
                wavefunct(:,:,ikx,iky,ikz)'* ...
                (dHdkz(:,:,ikx,iky,ikz)*sx + ...
                sx*dHdkz(:,:,ikx,iky,ikz)) * ...
                wavefunct(:,:,ikx,iky,ikz)/2;
            Asmnz(:,:,ikx,iky,ikz) = ...
                wavefunct(:,:,ikx,iky,ikz)'* ...
                (dHdkx(:,:,ikx,iky,ikz)*sy + ...
                sy*dHdkx(:,:,ikx,iky,ikz)) * ...
                wavefunct(:,:,ikx,iky,ikz)/2;
            for iband = 1:nw
            for jband = 1:nw
                endiff = energy(iband,ikx,iky,ikz) ...
                    -energy(jband,ikx,iky,ikz);
                if abs(endiff) > en_deg
                Amnx(iband,jband,ikx,iky,ikz) = ...
                   Amnx(iband,jband,ikx,iky,ikz)/endiff/1i;
                Amny(iband,jband,ikx,iky,ikz) = ...
                   Amny(iband,jband,ikx,iky,ikz)/endiff/1i;
                Amnz(iband,jband,ikx,iky,ikz) = ...
                   Amnz(iband,jband,ikx,iky,ikz)/endiff/1i;
               Asmnx(iband,jband,ikx,iky,ikz) = ...
                   Asmnx(iband,jband,ikx,iky,ikz)/endiff/1i;
               Asmny(iband,jband,ikx,iky,ikz) = ...
                   Asmny(iband,jband,ikx,iky,ikz)/endiff/1i;
               Asmnz(iband,jband,ikx,iky,ikz) = ...
                   Asmnz(iband,jband,ikx,iky,ikz)/endiff/1i;
                else
                    Amnx(iband,jband,ikx,iky,ikz) = 0;
                    Amny(iband,jband,ikx,iky,ikz) = 0;
                    Amnz(iband,jband,ikx,iky,ikz) = 0;
                    Asmnx(iband,jband,ikx,iky,ikz) = 0;
                    Asmny(iband,jband,ikx,iky,ikz) = 0;
                    Asmnz(iband,jband,ikx,iky,ikz) = 0;
                end
            end
            end
            
% Intraband (spin) Berry connection. Ann = i * < u_n | d/dk | u_n > ...
% = \sum_m < u_n | Amn | u_m >, Asnn = \sum_m < u_n | Asnn | u_m >. This
% part is not tested.
            for iband = 1:nw
                for jband = 1:nw
                    endiff = energy(iband,ikx,iky,ikz) ...
                    -energy(jband,ikx,iky,ikz);
                if abs(endiff) > en_deg
                    Amnx(iband,iband,ikx,iky,ikz) = ...
                        Amnx(iband,iband,ikx,iky,ikz) +...
                        wavefunct(:,iband,ikx,iky,ikz)' ...
                        *Amnx(:,:,ikx,iky,ikz) ...
                        *wavefunct(:,jband,ikx,iky,ikz);
                   Amny(iband,iband,ikx,iky,ikz) = ...
                        Amny(iband,iband,ikx,iky,ikz) +...
                        wavefunct(:,iband,ikx,iky,ikz)' ...
                        *Amny(:,:,ikx,iky,ikz) ...
                        *wavefunct(:,jband,ikx,iky,ikz);
                   Amnz(iband,iband,ikx,iky,ikz) = ...
                        Amnz(iband,iband,ikx,iky,ikz) +...
                        wavefunct(:,iband,ikx,iky,ikz)' ...
                        *Amnz(:,:,ikx,iky,ikz) ...
                        *wavefunct(:,jband,ikx,iky,ikz);
                   Asmnx(iband,iband,ikx,iky,ikz) = ...
                        Asmnx(iband,iband,ikx,iky,ikz) +...
                        wavefunct(:,iband,ikx,iky,ikz)' ...
                        *Asmnx(:,:,ikx,iky,ikz) ...
                        *wavefunct(:,jband,ikx,iky,ikz);
                   Asmny(iband,iband,ikx,iky,ikz) = ...
                        Asmny(iband,iband,ikx,iky,ikz) +...
                        wavefunct(:,iband,ikx,iky,ikz)' ...
                        *Asmny(:,:,ikx,iky,ikz) ...
                        *wavefunct(:,jband,ikx,iky,ikz);
                   Asmnz(iband,iband,ikx,iky,ikz) = ...
                        Asmnz(iband,iband,ikx,iky,ikz) +...
                        wavefunct(:,iband,ikx,iky,ikz)' ...
                        *Asmnz(:,:,ikx,iky,ikz) ...
                        *wavefunct(:,jband,ikx,iky,ikz);
                end
                end
            end
% End of intraband Berry connection

        end
    end
end

end


function [Amnx,Amny,Amnz,Asmnx,Asmny,Asmnz] = ...
    udu(wavefunct,kx,ky,kz,dkx,dky,dkz,blatt,nw)
% We calculate Berry connection using < u_n | d/dk | u_m >. This subroutine
% is not well-developed... Since there is a gauge problem in d/dk | u_m >,
% then a direct gradient is not appropriate. In order to solve this
% problem, one can use d/dk|u_n> = (e^-{i\theta}u_n,k+dk - u_n,k)/{dk},
% where e^{i\theta} = <u_n,k+dk | u_n,k> / |<u_n,k+dk | u_n,k>|.
gk1_un = zeros(nw,nw,length(kx),length(ky),length(kz));
gk2_un = zeros(nw,nw,length(kx),length(ky),length(kz));
gk3_un = zeros(nw,nw,length(kx),length(ky),length(kz));
gkx_un = zeros(nw,nw,length(kx),length(ky),length(kz));
gky_un = zeros(nw,nw,length(kx),length(ky),length(kz));
gkz_un = zeros(nw,nw,length(kx),length(ky),length(kz));
Amnx = zeros(nw,nw,length(kx),length(ky),length(kz));
Amny = zeros(nw,nw,length(kx),length(ky),length(kz));
Amnz = zeros(nw,nw,length(kx),length(ky),length(kz));
Asmnx = zeros(nw,nw,length(kx),length(ky),length(kz));
Asmny = zeros(nw,nw,length(kx),length(ky),length(kz));
Asmnz = zeros(nw,nw,length(kx),length(ky),length(kz));
for iband = 1:nw
    for ikx = 1:length(kx)
    for iky = 1:length(ky)
    for ikz = 1:length(kz)
        if ikx == 1
            e2itheta = wavefunct(:,iband,ikx+1,iky,ikz)' ...
                *wavefunct(:,iband,ikx,iky,ikz);
            e2itheta = e2itheta/abs(e2itheta);
            k1un = e2itheta^-1 * wavefunct(:,iband,ikx+1,iky,ikz) ...
                - wavefunct(:,iband,ikx,iky,ikz);
            k1un = k1un/dkx;
        elseif ikx == length(kx)
            e2itheta = wavefunct(:,iband,ikx,iky,ikz)' ...
                *wavefunct(:,iband,ikx-1,iky,ikz);
            e2itheta = e2itheta/abs(e2itheta);
            k1un = e2itheta^-1 * wavefunct(:,iband,ikx,iky,ikz) ...
                - wavefunct(:,iband,ikx-1,iky,ikz);
            k1un = k1un/dkx;
        else
            e2itheta1 = wavefunct(:,iband,ikx+1,iky,ikz)' ...
                *wavefunct(:,iband,ikx,iky,ikz);
            e2itheta1 = e2itheta1/abs(e2itheta1);
            e2itheta2 = wavefunct(:,iband,ikx,iky,ikz)' ...
                *wavefunct(:,iband,ikx-1,iky,ikz);
            e2itheta2 = e2itheta2/abs(e2itheta2);
            k1un = e2itheta1^-1 * wavefunct(:,iband,ikx+1,iky,ikz) ...
                - wavefunct(:,iband,ikx,iky,ikz) + ...
                e2itheta2^-1 * wavefunct(:,iband,ikx,iky,ikz) ...
                - wavefunct(:,iband,ikx-1,iky,ikz);
            k1un = k1un/dkx/2;
        end
        if iky == 1
            e2itheta = wavefunct(:,iband,ikx,iky+1,ikz)' ...
                *wavefunct(:,iband,ikx,iky,ikz);
            e2itheta = e2itheta/abs(e2itheta);
            k2un = e2itheta^-1 * wavefunct(:,iband,ikx,iky+1,ikz) ...
                - wavefunct(:,iband,ikx,iky,ikz);
            k2un = k2un/dky;
        elseif iky == length(ky)
            e2itheta = wavefunct(:,iband,ikx,iky,ikz)' ...
                *wavefunct(:,iband,ikx,iky-1,ikz);
            e2itheta = e2itheta/abs(e2itheta);
            k2un = e2itheta^-1 * wavefunct(:,iband,ikx,iky,ikz) ...
                - wavefunct(:,iband,ikx,iky-1,ikz);
            k2un = k2un/dky;
        else
            e2itheta1 = wavefunct(:,iband,ikx,iky+1,ikz)' ...
                *wavefunct(:,iband,ikx,iky,ikz);
            e2itheta1 = e2itheta1/abs(e2itheta1);
            e2itheta2 = wavefunct(:,iband,ikx,iky,ikz)' ...
                *wavefunct(:,iband,ikx,iky-1,ikz);
            e2itheta2 = e2itheta2/abs(e2itheta2);
            k2un = e2itheta1^-1 * wavefunct(:,iband,ikx,iky+1,ikz) ...
                - wavefunct(:,iband,ikx,iky,ikz) + ...
                e2itheta2^-1 * wavefunct(:,iband,ikx,iky,ikz) ...
                - wavefunct(:,iband,ikx,iky-1,ikz);
            k2un = k2un/dky/2;
        end
        if ikz == 1
            e2itheta = wavefunct(:,iband,ikx,iky,ikz+1)' ...
                *wavefunct(:,iband,ikx,iky,ikz);
            e2itheta = e2itheta/abs(e2itheta);
            k3un = e2itheta^-1 * wavefunct(:,iband,ikx,iky,ikz+1) ...
                - wavefunct(:,iband,ikx,iky,ikz);
            k3un = k3un/dkz;
        elseif ikz == length(kz)
            e2itheta = wavefunct(:,iband,ikx,iky,ikz)' ...
                *wavefunct(:,iband,ikx,iky,ikz-1);
            e2itheta = e2itheta/abs(e2itheta);
            k3un = e2itheta^-1 * wavefunct(:,iband,ikx,iky,ikz) ...
                - wavefunct(:,iband,ikx,iky,ikz-1);
            k3un = k3un/dkz;
        else
            e2itheta1 = wavefunct(:,iband,ikx,iky,ikz+1)' ...
                *wavefunct(:,iband,ikx,iky,ikz);
            e2itheta1 = e2itheta1/abs(e2itheta1);
            e2itheta2 = wavefunct(:,iband,ikx,iky,ikz)' ...
                *wavefunct(:,iband,ikx,iky,ikz-1);
            e2itheta2 = e2itheta2/abs(e2itheta2);
            k3un = e2itheta1^-1 * wavefunct(:,iband,ikx,iky,ikz+1) ...
                - wavefunct(:,iband,ikx,iky,ikz) + ...
                e2itheta2^-1 * wavefunct(:,iband,ikx,iky,ikz) ...
                - wavefunct(:,iband,ikx,iky,ikz-1);
            k3un = k3un/dkz/2;
        end
        
        gk1_un(:,iband,ikx,iky,ikz) = k1un;
        gk2_un(:,iband,ikx,iky,ikz) = k2un;
        gk3_un(:,iband,ikx,iky,ikz) = k3un;
    end
    end
    end
end
% Since the derivatives gk1_un ... gk3_un are along G_1, G_2, and G_3, we
% need to rotate them to xyz.
n1 = blatt(1,:);
n2 = blatt(2,:);
n3 = blatt(3,:);
n1 = n1/dot(n1,n1);
n2 = n2/dot(n2,n2);
n3 = n3/dot(n3,n3);
b = [n1;n2;n3];

for ikx = 1:length(kx)
for iky = 1:length(ky)
for ikz = 1:length(kz)
    for iband = 1:nw
    for jband = 1:nw
        A = b \ [gk1_un(iband,jband,ikx,iky,ikz);...
            gk2_un(iband,jband,ikx,iky,ikz); ...
            gk3_un(iband,jband,ikx,iky,ikz)];
        gkx_un(iband,jband,ikx,iky,ikz) = A(1);
        gky_un(iband,jband,ikx,iky,ikz) = A(2);
        gkz_un(iband,jband,ikx,iky,ikz) = A(3);
    end
    end
end
end
end

for ikx = 1:length(kx)
for iky = 1:length(ky)
for ikz = 1:length(kz)
    Amnx(:,:,ikx,iky,ikz) = ...
       1i*wavefunct(:,:,ikx,iky,ikz)' * gkx_un(:,:,ikx,iky,ikz);
   Amny(:,:,ikx,iky,ikz) = ...
       1i*wavefunct(:,:,ikx,iky,ikz)' * gky_un(:,:,ikx,iky,ikz);         
    Amnz(:,:,ikx,iky,ikz) = ...
       1i*wavefunct(:,:,ikx,iky,ikz)' * gkz_un(:,:,ikx,iky,ikz);
end
end
end

end


function [R,HR,wtR,wc,nw,nR] = getwanhr(head)
% One has to split wannier90_hr.dat file into two parts: One part is called
% head+wtR.dat, which only contains the fourth line part; The other part is
% called head+HmnR.dat which contains the [R1, R2, R3, mband, nband, Hr_r, Hr_i]
% data. A third file is head+wcenter.dat, which only contains wanner
% center coordinates
fileID = fopen(join([head,'wtR.dat']),'r'); % if head = '', then only wtR.dat
Rdeg = fscanf(fileID,'%f');
nR = length(Rdeg);
fclose(fileID);
wtR = 1./Rdeg;  % weight of each R vectors
fileID = fopen(join([head,'HmnR.dat']),'r');
C = textscan(fileID,'%d %d %d %d %d %n %n');
R1 = C{1};
R2 = C{2};
R3 = C{3};
Hr1 = C{6};
Hr2 = C{7};
nlines = length(R1);
nw = sqrt(nlines/nR);
fclose(fileID);
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

%fildID = fopen(join([head,'wcenter.dat']),'r');
%C = textscan(fildID,'%f %f %f');
%wc1 = C{1};
%wc2 = C{2};
%wc3 = C{3};
%fclose(fileID);
wc = zeros(nw,3);
%for i = 1:nw
%    wc(i,1) = wc1(i);
%    wc(i,2) = wc2(i);
%    wc(i,3) = wc3(i);
%end

end

function [Hk,dHdkx,dHdky,dHdkz] = buildHk(R,HR,wtR,wc,nw,nR,kpt,alatt)
% Note that here kpoint coordinate is in direct mode (in unit of reciprocal
% space
Hk = zeros(nw,nw);
dHdkx = zeros(nw,nw);
dHdky = zeros(nw,nw);
dHdkz = zeros(nw,nw);
for iR = 1:nR
    kdotR = 2*pi*R(iR,:)*kpt';
    aR = R(iR,:)*alatt;
    for i = 1:nw
    for j = 1:nw
    Hk(i,j) = Hk(i,j) + wtR(iR) * HR(iR,i,j) * exp(1i*kdotR);
    dHdkx(i,j) = dHdkx(i,j) + 1i*wtR(iR)*aR(1)*HR(iR,i,j)*exp(1i*kdotR);
    dHdky(i,j) = dHdky(i,j) + 1i*wtR(iR)*aR(2)*HR(iR,i,j)*exp(1i*kdotR);
    dHdkz(i,j) = dHdkz(i,j) + 1i*wtR(iR)*aR(3)*HR(iR,i,j)*exp(1i*kdotR);
    end
    end
end
Hk = roundn(Hk,-10);
dHdkx = roundn(dHdkx,-10);
dHdky = roundn(dHdky,-10);
dHdkz = roundn(dHdkz,-10);
end














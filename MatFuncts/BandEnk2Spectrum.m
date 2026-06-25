function SpectE = BandEnk2Spectrum(Enk,Znk)
% This function can simulate a typical band structure from Enk dataset to
% convert it in something like an APRES spectrum. The first input parameter
% is energy band, in the size of Enk(nkpt, nband). The second input
% parameter is optional, in the same size. It gives weight of each state.
% If not specified, then we just give a uniform value for each state.

E_exd = 0.5; % adding 500 meV above and below the lowest and highest band energy
delta_E = 0.001;  % energy step, in eV. Typical value is 1 meV
E_smear = 0.05;  % Smearing factor in the delta function, in eV
nk = size(Enk,1);
nband = size(Enk,2);
if nargin < 2
    Znk = ones(nk,nband);
end
Emin = min(min(Enk));
Emax = max(max(Enk));
Zmax = max(max(abs(Znk)));
if Zmax < 1e-3
    warning('The color scheme value might be too small')
end
Znk = Znk/Zmax; % To scale the color scheme value between [-1, 1]
ne = floor((Emax - Emin+2*E_exd)/delta_E);
SpectE = zeros(nk,ne);
for ik = 1:nk
    for ie = 1:ne
        e = Emin - E_exd + delta_E * (ie-1);
        for iband = 1:nband
        smear_fac = max(E_smear * abs(Znk(ik,iband)),0.01); % scales with color value
        % smear_fac = E_smear; % no scaling, with a uniform value
        SpectE(ik,ie) = SpectE(ik,ie) + Znk(ik,iband) ...
            * delta_funct(Enk(ik,iband)-e,smear_fac);
        end
    end
end

% plotting SpectE
X = linspace(1,nk,nk);
Y = linspace(Emin-E_exd, Emin-E_exd+ne*delta_E, ne);
[X, Y] = meshgrid(X,Y);
surface(X', Y', SpectE);
shading interp
colormap(flipud(cbrewer('seq','YlGnBu',256,'linear')))
colorbar
axis([1 nk Emin-E_exd Emin-E_exd+ne*delta_E])
ylabel('Energy (eV)','FontSize',28,'FontName','times new roman')
box on
ax = gca;
ax.LineWidth = 2;
ax.FontSize = 28;
ax.FontName = 'Times New Roman';
set(gca,'xticklabel',{[]})
xlabel('\it k-\rm path','FontSize',28,'FontName','times new roman')
% hold on
% hline = plot([1 nk],[0 0]); % Fermi level
% set(hline,'Color','w','LineStyle','--','LineWidth',1.5);

end
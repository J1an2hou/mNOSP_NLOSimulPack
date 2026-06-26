function Eg = plotting_banddispersion(kcorr,Enk,EF,nw,nkpts,KLABEL,klines,ksep,sz_k)
figure('Color', 'w');
hold on
if nargin < 9
for iw = 1:nw
    plot(kcorr', Enk(:,iw)-EF, 'LineWidth',2,'Color',[0 0 0])
end
else
for iw = 1:nw
    scatter(kcorr', Enk(:,iw)-EF,60,sz_k(:,iw),'filled','LineWidth',0.5,'MarkerEdgeColor','k')
    
end
colorbar
colormap(flipud(cbrewer('div','RdYlBu',256,'linear')))
end

ylabel('Energy (eV)','fontname','times new roman','fontsize',20)
hline = plot([0 kcorr(nkpts)],[0 0]); % Fermi level
set(hline,'Color','k','LineStyle','--')
kdash = zeros(1,klines+1);
for ik = 1:klines+1
    kdash(ik) = kcorr((ik-1)*ksep-ik+2);
end
xlim([kcorr(1), kcorr(end)])

yuplim = max(max(Enk)) - EF;
ydownlim = min(min(Enk)) - EF;
ycross = yuplim - ydownlim;
ymargin = 0.1*ycross;
yup = yuplim + ymargin;
ydown = ydownlim - ymargin;
ylim([ydown yup])
xticks(kdash)
xticklabels(KLABEL)
set(gca,'TickLabelInterpreter', 'latex','fontname','times new roman','fontsize',20)
set(gca,'linewidth',1.5)

box on
% Next, evaluate band gap
% Flatten the matrix to handle any dimension
E_flat = Enk(:) - EF;

% Find largest negative
negative_mask = E_flat < 0;
if any(negative_mask)
    largest_negative = max(E_flat(negative_mask));
else
    largest_negative = NaN;
    disp('Check Fermi level value\n');
end

% Find smallest positive
positive_mask = E_flat > 0;
if any(positive_mask)
    smallest_positive = min(E_flat(positive_mask));
else
    smallest_positive = NaN;
    fprintf('Check Fermi level value\n');
end

Eg = smallest_positive - largest_negative;
disp(['Band gap is ',num2str(Eg)])
end
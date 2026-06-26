function plt_surf = plotting_surf_onXY(xx,yy,nx,ny,Zk)
% the xx and yy and Zk are column vectors with same sizes
% nx and ny are required grid point that can be very large
figure('Color', 'w')
X = linspace(min(xx),max(xx),nx);
Y = linspace(min(yy),max(yy),ny);
[X,Y] = meshgrid(X,Y);
Z = griddata(xx',yy',Zk,X,Y,'linear');
surf(X, Y, Z)
shading interp
colorbar
colormap(cbrewer('div','RdYlBu',512,'linear'))
xlabel('$k_x$','interpreter','latex','fontsize',20,'FontName','times new roman')
ylabel('$k_y$','interpreter','latex','fontsize',20,'FontName','times new roman')
box on
set(gca,'FontSize',20)
set(gca,'LineWidth',1.5)
set(gca,'FontName','Times New Roman')
set(gca,'Color','w')
plt_surf = true;
% axis equal
axis tight
end
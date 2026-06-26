function plt_quiver = plotting_quiver_onXY(xx, yy, nx, ny, Z1, Z2, Z3, cmap_name)
% the xx and yy are column vectors with same sizes
% Z1 and Z2 are the two components (e.g., x and y components of vectors)
% Z3 is an optional parameter for background color gradient
% cmap_name is optional colormap name (default: 'RdYlBu')
% nx and ny are required grid points that can be very large

% Set default colormap if not provided
if nargin < 8 || isempty(cmap_name)
    cmap_name = 'RdYlBu';
end

figure('Color', 'w')

% Create grid
X = linspace(min(xx), max(xx), nx);
Y = linspace(min(yy), max(yy), ny);
[X, Y] = meshgrid(X, Y);

% Interpolate vector components onto the grid
U = griddata(xx', yy', Z1, X, Y, 'linear');  % X-component of vectors
V = griddata(xx', yy', Z2, X, Y, 'linear');  % Y-component of vectors

% Check if Z3 is provided
if nargin >= 7 && ~isempty(Z3)
    % Interpolate Z3 for background
    Z_bg = griddata(xx', yy', Z3, X, Y, 'linear');
    
    % Plot background first (so quiver appears on top)
    surf(X, Y, zeros(size(X)), Z_bg, 'EdgeColor', 'none')
    view(2)  % 2D view from top
    shading interp
    colorbar
    
    % Apply colormap
    try
        colormap(cbrewer('div', cmap_name, 512, 'linear'))
    catch
        % Fallback to default MATLAB colormap if cbrewer fails
        colormap(jet)
    end
    
    hold on
    
    % Add quiver plot on top
    q = quiver(X, Y, U, V, 'k', 'LineWidth', 1.5, 'AutoScale', 'on');
    
    % Optional: Make quiver arrows more visible
    q.MaxHeadSize = 0.3;  % Adjust arrow head size
    q.AutoScaleFactor = 0.8;  % Adjust arrow scaling
else
    % No background - just quiver plot
    q = quiver(X, Y, U, V, 'k', 'LineWidth', 1.5, 'AutoScale', 'on');
    % Optional: Make quiver arrows more visible
    q.MaxHeadSize = 0.3;  % Adjust arrow head size
    q.AutoScaleFactor = 0.8;  % Adjust arrow scaling
    % colorbar  % Remove if you don't want colorbar
end

xlabel('$k_x$', 'interpreter', 'latex', 'fontsize', 20, 'FontName', 'times new roman')
ylabel('$k_y$', 'interpreter', 'latex', 'fontsize', 20, 'FontName', 'times new roman')
box on
set(gca, 'FontSize', 20)
set(gca, 'LineWidth',1.5)
set(gca, 'FontName', 'Times New Roman')
set(gca, 'Color', 'w')
axis equal
axis tight

plt_quiver = true;
end
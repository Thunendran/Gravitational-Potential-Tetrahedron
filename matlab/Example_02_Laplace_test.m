% =========================================================================
% Example 02: Laplacian Validation for the Gravitational Potential
% -------------------------------------------------------------------------
% Description:
%     This example verifies the Laplace equation for the gravitational
%     potential of a homogeneous tetrahedron using the analytic formulation.
%
%     Specifically, it numerically evaluates the Laplacian (∇²V) of the
%     potential in the plane z = 1/2, confirming:
%
%           ∇²V = -4π   → Inside the tetrahedron (Poisson equation)
%           ∇²V =  0    → Outside the tetrahedron (Laplace equation)
%
%     The result is presented as a 2D color map of the Laplacian values,
%     with the tetrahedral boundary overlaid as a white polygon.
%
% Reference:
%     Periyandy, T. & Bevis, M. (2025)
%     “The Gravitational Potential Inside, On and Outside of a Homogeneous
%      Tetrahedron.”
%
% Authors:
%     Thunendran Periyandy  (corresponding author)
%     Michael Bevis
%
% Date:    July 31, 2025
% Version: 1.0
% =========================================================================

close all
clear all
clc

% -------------------------------------------------------------------------
% Define Tetrahedron Vertices
% -------------------------------------------------------------------------
% The tetrahedron is defined as the unit simplex with one vertex at the
% origin and edges along the coordinate axes.
%
% Vertices:
%   A = (0, 0, 0)
%   B = (0, 1, 0)
%   C = (1, 0, 0)
%   D = (0, 0, 1)
% -------------------------------------------------------------------------
vertices = [ 0.0 0.0 0.0;
             0.0 1.0 0.0;
             1.0 0.0 0.0;
             0.0 0.0 1.0 ];
vA = vertices(1,:); vB = vertices(2,:); vC = vertices(3,:); vD = vertices(4,:);

% -------------------------------------------------------------------------
% Create Laplace Tester Object
% -------------------------------------------------------------------------
% The `TetrahedronLaplaceTest` class numerically evaluates the Laplacian of
% the analytical potential on a specified grid.
% -------------------------------------------------------------------------
laplace_tester = TetrahedronLaplaceTest(vA, vB, vC, vD);

% -------------------------------------------------------------------------
% Compute Laplacian Grid
% -------------------------------------------------------------------------
% Computes ∇²V numerically across a uniform grid in the plane z = 1/2.
%
% The grid resolution (here, 200 × 200) controls the smoothness of the
% visualization.
% -------------------------------------------------------------------------
z_fixed = 0.5;
x_range = [-0.1, 0.6];
y_range = [-0.1, 0.6];
[X,Y,laplacian_grid] = laplace_tester.compute_laplacian_grid(z_fixed, x_range, y_range, 200);

% -------------------------------------------------------------------------
% Compute Edge–Plane Intersections
% -------------------------------------------------------------------------
% Determines where tetrahedral edges intersect the plane z = 1/2.
% These intersection points define the polygonal boundary of the
% tetrahedron’s cross-section for plotting.
% -------------------------------------------------------------------------
edge_intersection = @(v1,v2,z_plane) ...
    ( (v1(3)-z_plane)*(v2(3)-z_plane) < 0 ) .* ...
    ( v1 + ((z_plane - v1(3)) / (v2(3) - v1(3))) * (v2 - v1) );

edges = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
intersection_points = [];
for e = 1:size(edges,1)
    pt = edge_intersection(vertices(edges(e,1),:), vertices(edges(e,2),:), z_fixed);
    if any(pt) && ~all(isnan(pt))
        intersection_points = [intersection_points; pt];
    end
end

% -------------------------------------------------------------------------
% Define Custom Colormap
% -------------------------------------------------------------------------
% Matches the colormap used in the Python visualization:
%   Blue → Cyan → Yellow → Red
%
% The Laplacian range is set from -4π (interior) to 0 (exterior).
% -------------------------------------------------------------------------
N = 50;
colors = [0.1 0.4 1]; % Blue (start)
for i = 2:N-1
    ratio = (i-2)/(N-3);
    colors = [colors; ratio 1 1-ratio];
end
colors = [colors; 1 0 0]; % Red (end)
vmin = -4*pi; vmax = 0;

% -------------------------------------------------------------------------
% Plot Laplacian Field
% -------------------------------------------------------------------------
fig = figure('Position',[100 100 800 600]);
imagesc(x_range, y_range, laplacian_grid);
set(gca,'YDir','normal');
caxis([vmin vmax]);
colormap(colors);

% -------------------------------------------------------------------------
% Add Colorbar
% -------------------------------------------------------------------------
cb = colorbar('Ticks',[-12 -10 -8 -6 -4 -2 0], ...
              'TickLabels',arrayfun(@num2str,[-12 -10 -8 -6 -4 -2 0],'UniformOutput',false));
cb.Label.String = '$\nabla^2 V$';
cb.Label.Interpreter = 'latex';

hold on;

% -------------------------------------------------------------------------
% Overlay Tetrahedral Cross-Section
% -------------------------------------------------------------------------
% Draws the polygonal boundary where the tetrahedron intersects z = 1/2.
% Annotates regions corresponding to ∇²V = -4π (interior) and ∇²V = 0 (exterior).
% -------------------------------------------------------------------------
if size(intersection_points,1) >= 3
    K = convhull(intersection_points(:,1), intersection_points(:,2));
    plot(intersection_points(K,1), intersection_points(K,2), ...
         'w-', 'LineWidth',2, 'DisplayName','Tetrahedral Boundary');
    
    interior_x = mean(intersection_points(:,1));
    interior_y = mean(intersection_points(:,2));
    text(interior_x, interior_y, '$\nabla^2 V = -4\pi$', ...
        'Color','w', 'FontSize',15, 'FontWeight','bold', ...
        'HorizontalAlignment','center', 'Interpreter','latex');
end

% Exterior annotation
outside_x = x_range(2) - 0.3;
outside_y = y_range(2) - 0.3;
text(outside_x, outside_y, '$\nabla^2 V = 0$', ...
    'Color','w', 'FontSize',24, 'FontWeight','bold', ...
    'HorizontalAlignment','center', 'Interpreter','latex');

% -------------------------------------------------------------------------
% Final Formatting
% -------------------------------------------------------------------------
title('Laplacian Test in $z = \frac{1}{2}$ Plane of the Tetrahedron', ...
    'FontSize',26, 'Interpreter','latex');
xlabel('$x$', 'FontSize',24, 'Interpreter','latex');
ylabel('$y$', 'FontSize',24, 'Interpreter','latex');

legend('Location','SouthEast','Color','none','TextColor','w');
grid off;

% -------------------------------------------------------------------------
% Export Figure
% -------------------------------------------------------------------------
% Saves the Laplacian test result as a high-resolution PNG (300 dpi)
% -------------------------------------------------------------------------
print(fig, 'Laplacian_Test_Matlab.png', '-dpng', '-r300');

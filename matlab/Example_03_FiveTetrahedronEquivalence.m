% =========================================================================
% Example 03: Equivalence of the Gravitational Acceleration of a Cuboid
%             and the Superposition of Five Tetrahedral Accelerations
% -------------------------------------------------------------------------
% Description:
%     This example reproduces the equivalence relationship between the
%     gravitational field of a cube and the superposition of five tetrahedral
%     fields, as described in Appendix B.
%
%     The gravitational acceleration components are computed numerically
%     for each tetrahedron using the OOP-based analytical potential solver
%     (`TetrahedronPotentialCalculatorOOP`), and the results are displayed
%     using quiver plots for both the individual and combined fields.
%
% Reference:
%     Periyandy, T. & Bevis, M. (2025)
%     “The Gravitational Potential Inside, On and Outside of a
%      Homogeneous Tetrahedron.”
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
% Define the 5 tetrahedra that compose the cube
% Each tetrahedron occupies one of the five unique orientations whose sum
% reproduces the gravitational field of the full cube.
% -------------------------------------------------------------------------
tetrahedrons = {
    [-1 -1  1;  1 -1 -1; -1  1 -1; -1 -1 -1];
    [-1 -1  1;  1  1  1; -1  1 -1; -1  1  1];
    [ 1 -1 -1;  1  1 -1; -1  1 -1;  1  1  1];
    [-1 -1  1;  1  1  1;  1 -1 -1;  1 -1  1];
    [-1 -1  1;  1 -1 -1; -1  1 -1;  1  1  1]
};

% -------------------------------------------------------------------------
% Define computational grid in the z = 0 plane
% -------------------------------------------------------------------------
x_vals = linspace(-2, 2, 25);
y_vals = linspace(-2, 2, 25);
[X, Y] = meshgrid(x_vals, y_vals);
z_fixed = 0;

% -------------------------------------------------------------------------
% Function to compute gravitational acceleration field for one tetrahedron
% -------------------------------------------------------------------------
compute_field = @(verts) ...
    local_compute_field(X, Y, z_fixed, verts);

% -------------------------------------------------------------------------
% Compute acceleration field for each tetrahedron
% -------------------------------------------------------------------------
num_tetra = numel(tetrahedrons);
fields = cell(num_tetra, 1);
polygons = cell(num_tetra, 1);

disp('Computing gravitational acceleration fields for 5 tetrahedra...');
for k = 1:num_tetra
    vertices = tetrahedrons{k};
    [U, V] = compute_field(vertices);
    fields{k} = {U, V};

    % Compute intersection polygon (cross-section at z = 0)
    polygons{k} = intersection_polygon(vertices, z_fixed);
end

% -------------------------------------------------------------------------
% Compute combined field by superposition
% -------------------------------------------------------------------------
U_total = zeros(size(X));
V_total = zeros(size(Y));
for k = 1:num_tetra
    U_total = U_total + fields{k}{1};
    V_total = V_total + fields{k}{2};
end

% -------------------------------------------------------------------------
% Create 2×3 subplot grid (5 tetrahedra + combined field)
% -------------------------------------------------------------------------
figure('Position', [100 100 1200 800]);
titles = arrayfun(@(i) sprintf('Tetrahedron %d', i), 1:5, 'UniformOutput', false);
titles{6} = 'Combined Field';

for idx = 1:6
    subplot(2, 3, idx);
    
    if idx <= 5
        quiver(X, Y, fields{idx}{1}, fields{idx}{2}, 'r', 'LineWidth', 1);
        hold on;
        poly_pts = polygons{idx};
    else
        quiver(X, Y, U_total, V_total, 'r', 'LineWidth', 1);
        hold on;
        poly_pts = vertcat(polygons{:});
    end

    % Draw boundary polygon
    if size(poly_pts, 1) >= 3
        K = convhull(poly_pts(:,1), poly_pts(:,2));
        plot(poly_pts(K,1), poly_pts(K,2), 'k-', 'LineWidth', 1.5);
    end
    
    title(titles{idx}, 'FontSize', 14, 'FontWeight', 'bold');
    xlim([-2 2]); ylim([-2 2]);
    axis equal; grid off;
end

sgtitle(['Equivalence of the Gravitational Acceleration of a Cuboid ' ...
         'and the Superposition of Five Tetrahedral Accelerations ' ...
         'in the Plane z = 0'], ...
         'FontSize', 16, 'FontWeight', 'bold');

%print('Five_Tetrahedrons_Field_Equivalence_Matlab.png', '-dpng', '-r400');
%disp('Field plots saved successfully.');

% =========================================================================
% --- Local Helper Function: Compute Field for a Single Tetrahedron ---
% =========================================================================
function [U, V] = local_compute_field(X, Y, z_fixed, vertices)
    % local_compute_field
    % Computes gravitational acceleration components (U,V) in the z = z_fixed
    % plane for the specified tetrahedron using the OOP solver.
    
    U = zeros(size(X));
    V = zeros(size(Y));
    
    vA = vertices(1,:); vB = vertices(2,:); vC = vertices(3,:); vD = vertices(4,:);
    
    for i = 1:size(X,1)
        for j = 1:size(X,2)
            P = [X(i,j), Y(i,j), z_fixed];
            calc = TetrahedronPotentialCalculatorOOP(vA, vB, vC, vD, P);
            
            % Numerical gradient for gravitational acceleration
            h = 1e-4;
            gx = -(calc.compute([P(1)+h,P(2),P(3)]) - calc.compute([P(1)-h,P(2),P(3)])) / (2*h);
            gy = -(calc.compute([P(1),P(2)+h,P(3)]) - calc.compute([P(1),P(2)-h,P(3)])) / (2*h);
            
            U(i,j) = gx;
            V(i,j) = gy;
        end
    end
end

% =========================================================================
% --- Local Helper Function: Compute Intersection Polygon ---
% =========================================================================
function poly_points = intersection_polygon(vertices, z_plane)
    % intersection_polygon
    % Computes the polygonal cross-section formed by the intersection of
    % the tetrahedron with the specified z-plane.
    
    edges = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
    poly_points = [];
    for e = 1:size(edges,1)
        v1 = vertices(edges(e,1),:);
        v2 = vertices(edges(e,2),:);
        if (v1(3)-z_plane)*(v2(3)-z_plane) < 0
            t = (z_plane - v1(3)) / (v2(3) - v1(3));
            pt = v1 + t*(v2 - v1);
            poly_points = [poly_points; pt];
        end
    end
end

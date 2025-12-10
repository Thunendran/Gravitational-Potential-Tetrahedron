%% Compare GPolyhedron vs TetrahedronPotentialCalculatorOOP for key points
clear; clc;

% Define 5 test points
points = [
    0, 0, 0;
    1, 0, 0;
    0, 1, 0;
    0, 0, 1;
    1/4, 1/4, 1/4
];
nP = size(points, 1);

% --- Constants ---
G = 1.0; rho = 1.0;

% --- Paths ---
projectRoot = fileparts(mfilename('fullpath'));
coordsFile = fullfile(projectRoot, 'GPolyhedron-main', 'coords.txt');
topoFile   = fullfile(projectRoot, 'GPolyhedron-main', 'topology.txt');

% --- Load tetrahedron vertices ---
A = readmatrix(fullfile(projectRoot, '..', 'data', 'vertices.txt'));
A1 = A(1,:); B = A(2,:); C = A(3,:); D = A(4,:);
calc = TetrahedronPotentialCalculatorOOP(A1,B,C,D,[],G,rho);

% --- Initialize results ---
V_tet = zeros(nP,1);
V_poly = zeros(nP,1);

% --- Loop through the five test points ---
for i = 1:nP
    P = points(i,:);
    % Compute potential using Periyandy @ Bevis code
    V_tet(i) = calc.compute(P);

    % Prepare shifted coords for GPolyhedron
    shifted = A - P;
    tmpCoords = fullfile(projectRoot, 'GPolyhedron-main', 'coords_tmp.txt');
    writematrix(shifted, tmpCoords, 'Delimiter', ' ');

    % Compute potential using Tsoulis (2021)
    [Vtmp, ~, ~, ~, ~] = GPolyhedron(G, rho, topoFile, tmpCoords);
    V_poly(i) = Vtmp;
end

% --- Build comparison table ---
tableData = [(1:nP)', points, V_tet, V_poly];
headers = {'ID','x','y','z','Potential (Periyandy & Bevis)','Potential (Tsoulis, 2021)'};

% --- Save as text file ---
fileID = fopen(fullfile(projectRoot, 'comparison_table.txt'), 'w');
fprintf(fileID, '%-4s %-10s %-10s %-10s %-30s %-30s\n', headers{:});
for i = 1:nP
    fprintf(fileID, '%-4d %-10.6f %-10.6f %-10.6f %-30.15e %-30.15e\n', ...
        i, points(i,1), points(i,2), points(i,3), V_tet(i), V_poly(i));
end
fclose(fileID);

fprintf('Comparison table saved as comparison_table.txt\n');

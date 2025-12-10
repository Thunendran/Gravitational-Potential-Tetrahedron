%% Compare GPolyhedron vs TetrahedronPotentialCalculatorOOP for key points + epsilon shifts
clear; clc;

% --- Small perturbation for stability test ---
epsilon = 1e-10;

% --- Define base test points ---
basePoints = [
    0, 0, 0;
    1, 0, 0;
    0, 1, 0;
    0, 0, 1;
    1/4, 1/4, 1/4
];
nBase = size(basePoints, 1);

% --- Generate epsilon-shifted points ---
% For each base point, add one perturbed version with +epsilon added to all coordinates
epsPoints = basePoints + epsilon;
points = [basePoints; epsPoints];
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

% --- Loop through all test points (original + perturbed) ---
for i = 1:nP
    P = points(i,:);
    % Compute potential using Periyandy & Bevis method
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
pointType = [repmat({'Base'}, nBase, 1); repmat({'+Eps'}, nBase, 1)];
tableData = [(1:nP)', points, V_tet, V_poly];
headers = {'ID','Type','x','y','z','Potential (Periyandy & Bevis)','Potential (Tsoulis, 2021)'};

% --- Save as text file ---
outFile = fullfile(projectRoot, 'comparison_table.txt');
fileID = fopen(outFile, 'w');
fprintf(fileID, '%-4s %-6s %-12s %-12s %-12s %-30s %-30s\n', headers{:});
for i = 1:nP
    fprintf(fileID, '%-4d %-6s %-12.10f %-12.10f %-12.10f %-30.15e %-30.15e\n', ...
        i, pointType{i}, points(i,1), points(i,2), points(i,3), V_tet(i), V_poly(i));
end
fclose(fileID);

fprintf('Comparison table with epsilon perturbations saved as %s\n', outFile);

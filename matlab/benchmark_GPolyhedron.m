%% Benchmark Test for GPolyhedron (Potential Only)
clear; clc;

% Constants
G = 1.0;   % gravitational constant
rho = 1.0; % density

% Paths
projectRoot = fileparts(mfilename('fullpath'));
coordsFile = fullfile(projectRoot, 'GPolyhedron-main', 'coords.txt');
topoFile   = fullfile(projectRoot, 'GPolyhedron-main', 'topology.txt');
pointsFile = fullfile(projectRoot, '..', 'data', 'points.txt');

% Check files
assert(isfile(coordsFile), 'coords.txt not found.');
assert(isfile(topoFile), 'topology.txt not found.');
assert(isfile(pointsFile), 'points.txt not found.');

% Load evaluation points
P = readmatrix(pointsFile);
nP = size(P, 1);

fprintf('Running GPolyhedron for %d points...\n', nP);

V = zeros(nP, 1);

tic;
for i = 1:nP
    % Shift coordinates so current point is origin
    A = readmatrix(coordsFile);
    shiftedCoords = A - P(i,:);
    tmpCoords = fullfile(projectRoot, 'GPolyhedron-main', 'coords_tmp.txt');
    writematrix(shiftedCoords, tmpCoords, 'Delimiter', ' ');
    
    % Compute potential only
    [Vtmp, ~, ~, ~, ~] = GPolyhedron(G, rho, topoFile, tmpCoords);
    V(i) = Vtmp;
end
elapsed = toc;

fprintf('GPolyhedron runtime: %.6f s\n', elapsed);

% Save results
writematrix(V, fullfile(projectRoot, 'GPolyhedron_results.txt'));
fid = fopen(fullfile(projectRoot, 'GPolyhedron_time.txt'), 'w');
fprintf(fid, '%.6e\n', elapsed);
fclose(fid);

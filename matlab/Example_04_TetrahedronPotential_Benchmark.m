% =========================================================================
% Example 04: Benchmarking the Analytical Gravitational Potential Solver
% -------------------------------------------------------------------------
% Description:
%     This script benchmarks the MATLAB implementation of the analytical
%     gravitational potential solver for a homogeneous tetrahedron using
%     the object-oriented class `TetrahedronPotentialCalculatorOOP`.
%
%     It reads the tetrahedron geometry and evaluation points from text
%     files, computes the potential values at all field points, measures
%     the total runtime, and writes the results to output files for
%     comparison against Julia and Python implementations.
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
% Date:    November 1, 2025
% Version: 1.0
% =========================================================================

clear; clc;

% -------------------------------------------------------------------------
% Step 1. Load Tetrahedron Geometry and Evaluation Points
% -------------------------------------------------------------------------
% The geometry file ('vertices.txt') contains four 3D vertex coordinates
% defining the tetrahedron, arranged as:
%
%       A1x  A1y  A1z
%       Bx   By   Bz
%       Cx   Cy   Cz
%       Dx   Dy   Dz
%
% The points file ('points.txt') contains an N×3 matrix of field points
% where the potential will be evaluated.
% -------------------------------------------------------------------------
A = readmatrix('../data/vertices.txt');   % 4×3 matrix of vertex coordinates
P = readmatrix('../data/points.txt');     % N×3 matrix of evaluation points

% -------------------------------------------------------------------------
% Step 2. Assign Vertex Labels for Clarity
% -------------------------------------------------------------------------
A1 = A(1,:); 
B  = A(2,:); 
C  = A(3,:); 
D  = A(4,:);

% -------------------------------------------------------------------------
% Step 3. Initialize Potential Calculator (OOP class)
% -------------------------------------------------------------------------
% Create an instance of the TetrahedronPotentialCalculatorOOP class
% The constructor automatically ensures correct orientation of vertices.
%
% Arguments:
%   A1, B, C, D : Vertex coordinates
%   []           : No fixed evaluation point (vectorized mode)
%   1.0, 1.0     : Gravitational constant G = 1, surface density σ = 1
% -------------------------------------------------------------------------
calc = TetrahedronPotentialCalculatorOOP(A1, B, C, D, [], 1.0, 1.0);

% -------------------------------------------------------------------------
% Step 4. Compute Potentials at All Points (Vectorized)
% -------------------------------------------------------------------------
% The compute_vectorized() method internally parallelizes the evaluation
% using MATLAB’s parallel pool (parfor), enabling efficient computation
% for large datasets (e.g., 100k points).
% -------------------------------------------------------------------------
tic;
U = calc.compute_vectorized(P);
elapsed = toc;

fprintf('MATLAB runtime: %.6f s\n', elapsed);

% -------------------------------------------------------------------------
% Step 5. Write Results to Output Files
% -------------------------------------------------------------------------
% The potential values are written to 'matlab_results.txt' (N×1 vector),
% and the total runtime is written to 'matlab_time.txt' for comparison
% with Julia and Python timings.
% -------------------------------------------------------------------------
writematrix(U, 'matlab_results.txt');

fid = fopen('matlab_time.txt', 'w');
fprintf(fid, '%.6e\n', elapsed);
fclose(fid);

fprintf('Results saved: matlab_results.txt\n');
fprintf('Runtime saved: matlab_time.txt\n');
fprintf('Benchmark completed successfully.\n');

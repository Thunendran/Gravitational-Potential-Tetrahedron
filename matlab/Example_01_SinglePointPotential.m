% =========================================================================
% Example 01: Gravitational Potential at a Single Evaluation Point
% -------------------------------------------------------------------------
% Description:
%     Demonstrates the computation of the gravitational potential generated
%     by a homogeneous tetrahedron at a single field point using the
%     object-oriented implementation `TetrahedronPotentialCalculatorOOP`.
%
%     This basic example corresponds to the analytical test case in which
%     the tetrahedron occupies the positive octant of a unit cube.
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
% Date:    July 31 2025
% Version: 1.0
% =========================================================================

% -------------------------------------------------------------------------
% Define tetrahedron vertices (A, B, C, D)
% The tetrahedron is a right-angled form with one vertex at (1, 1, 1)
% and edges aligned with coordinate axes.
% -------------------------------------------------------------------------
vA = [1.0, 1.0, 1.0];
vB = [2.0, 1.0, 1.0];
vC = [1.0, 2.0, 1.0];
vD = [1.0, 1.0, 2.0];

% -------------------------------------------------------------------------
% Define the evaluation point
% Here, the potential is evaluated at the origin (0, 0, 0)
% -------------------------------------------------------------------------
p_eval = [0.0, 0.0, 0.0];

% -------------------------------------------------------------------------
% Create a tetrahedron potential calculator object
% Inputs: vertices (A–D) and the evaluation point
% -------------------------------------------------------------------------
calc = TetrahedronPotentialCalculatorOOP(vA, vB, vC, vD, p_eval);

% -------------------------------------------------------------------------
% Compute gravitational potential at the evaluation point
% -------------------------------------------------------------------------
potential = calc.compute();

% -------------------------------------------------------------------------
% Display the numerical result with high precision
% -------------------------------------------------------------------------
fprintf('Gravitational Potential at point (0, 0, 0): %.15f\n', potential);


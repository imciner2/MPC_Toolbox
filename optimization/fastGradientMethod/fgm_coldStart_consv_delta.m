function [ delta ] = fgm_coldStart_consv_delta( L, N, G, g )
%FGM_COLDSTART_CONSV_DELTA Compute the conservative value for delta assuming cold-starting
%
% This function will compute the upper-bound for the Delta value used to
% create the lower iteration bound for the Fast Gradient Method.
%
% This bound is derived in: 
%   S. Richter, C. N. Jones, and M. Morari, “Computational Complexity
%   Certification for Real-Time MPC With Input Constraints Based on the
%   Fast Gradient Method,” IEEE Transactions on Automatic Control,
%   vol. 57, no. 6, pp. 1391–1403, 2012.
%
%
% Usage:
%   [ delta ] = FGM_COLDSTART_CONSV_DELTA( L, N, G, g);
%
% Inputs:
%   L    - The largest eigenvalue of the Hessian matrix
%   N    - The horizon length
%   G    - The input constraint matrix
%   g    - The input constraint vector
%   H    - The Hessian matrix
%   M    - The linear term from the cost function
%
% Output:
%   delta - The delta parameter
%
%
% Created by: Ian McInerney
% Created on: August 13, 2018
% Version: 1.0
% Last Modified: August 13, 2018
%
% Revision History
%   1.0 - Initial release


%% Find some variable sizing information
[n, m] = size(G);
I = eye(N*m);
Z = zeros(N*m, 1);


%% Create the options for quadprog
opt = optimoptions('quadprog');
opt.Display = 'off';


%% Find the vector that has the largest 2-norm in the admissible input set
Gl = kron( eye(N), G);
gl = kron( ones(N,1), g);

r = quadprog(I, Z, Gl, gl, [], [], [], [], [], opt);
r = norm(r);


%% Compute the conservative value of delta
delta = L/2 * r^2;

end
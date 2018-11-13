function [ G ] = condensed_primal_constraint_coefficient_gen( N, A, B, E, varargin )
%CONDENSED_PRIMAL_CONSTRAINT_COEFFICIENT_GEN Generate the constraint matrix G
%
% Create the condensed variable coefficient matrix for the inequality 
% constraints of the condensed linear time-invariant MPC problem.
%
%
% Usage:
%   [ G ] = CONDENSED_PRIMAL_CONSTRAINT_COEFFICIENT_GEN( N, A, B, E )
%   [ G ] = CONDENSED_PRIMAL_CONSTRAINT_COEFFICIENT_GEN( N, A, B, E, D )
%
% Inputs:
%   N - The horizon length
%   A - The state transition matrix
%   B - The input mapping matrix
%   E - The stage constraints for the inputs
%   D - The stage constraints for the states
%
% Outputs:
%   G - The coefficient matrix
%
% See also 
%
% Created by: Ian McInerney
% Created on: August 17, 2018
% Version: 1.0
% Last Modified: August 17, 2018
%
% Revision History
%   1.0 - Initial release  


%% Parse the input arguments
p = inputParser;
addOptional(p, 'D', []);
parse(p,varargin{:});

% Extract the matrices
D = p.Results.D;


%% Find the number of constraints and the system size
[nE, ~] = size(E);
[nD, ~] = size(D);
[n, m] = size(B);


%% Get the condensed prediction matrix
Gamma = condensed_prediction_gen(A, B, N);


%% Create the component matrices
Dbar = kron( eye(N), [D; zeros(nE, n)]);
Ebar = kron( eye(N), [zeros(nD,m); E]);


%% Put it all together
G = Dbar*Gamma + Ebar;

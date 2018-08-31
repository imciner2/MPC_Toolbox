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


[n, m] = size(B);

%% Parse the input arguments
p = inputParser;
addOptional(p, 'D', []);
parse(p,varargin{:});

% Extract the matrices
D = p.Results.D;


%% Get the condensed system matrices
Gamma = condensed_prediction_gen(A, B, N);
Phi = condensed_initial_gen(A, N);

% Shift the Gamma matrix
Gamma_shift = [eye(n*(N-1)), zeros(n*(N-1), n)]*Gamma;
Gamma_shift = [zeros(n, N*m);
               Gamma_shift];


%% Create the matrices
if ( isempty(D) )
    % If no state constraints are present
    G = kron(eye(N), E);
    return
elseif ( isempty(E) )
    % No input constraints are present
    Ebar = kron( eye(N), D*B );
    Dbar = kron( eye(N), D*A );
else
    [nE, ~] = size(E);
    Ebar = kron( eye(N), [D*B; E] );
    Dbar = kron( eye(N), [D*A; zeros(nE, n)] );
end


%% Put it all together
G = Dbar*Gamma_shift + Ebar;

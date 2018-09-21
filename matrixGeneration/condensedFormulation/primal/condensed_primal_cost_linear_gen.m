function [ J ] = condensed_primal_cost_linear_gen( N, A, B, Q, varargin )
%CONDENSED_PRIMAL_COST_LINEAR_GEN Generate the linear term matrix J
%
% Create the matrix for the linear term in the cost function of the condensed 
% linear time-invariant MPC problem.
%
% If no P matrix is supplied, it defaults to Q.
% If no S matrix is supplied, it defaults to 0.
%
% Usage:
%   [ J ] = CONDENSED_PRIMAL_COST_LINEAR_GEN( N, A, B, Q );
%   [ J ] = CONDENSED_PRIMAL_COST_LINEAR_GEN( N, A, B, Q, P );
%   [ J ] = CONDENSED_PRIMAL_COST_LINEAR_GEN( N, A, B, Q, P, S );
%
% Inputs:
%   N - The horizon length
%   A - The state transition matrix
%   B - The input mapping matrix
%   Q - The state weighting matrix
%   P - The final state weighting matrix
%   S - The state-input cross term weight matrix
%
% Outputs:
%   G - The linear term matrix
%
% See also CONDENSED_PRIMAL_COST_HESSIAN_GEN
%
% Created by: Ian McInerney
% Created on: May 22, 2018
% Version: 1.1
% Last Modified: August 17, 2018
%
% Revision History
%   1.0 - Initial release  
%   1.1 - Update the name, inputs & documentation
%   1.2 - Added S term


%% Parse the input arguments
p = inputParser;
addOptional(p, 'P', Q);
addOptional(p, 'S', zeros(n,m));
parse(p,varargin{:});

% Extract the matrices
P = p.Results.P;
S = p.Results.S;


%% Create the diagonal matrix for the Q matrix
Qbar = kron(eye(N-1), Q);
Qbar = blkdiag(Qbar, P);

Sbar = kron(eye(N), S);


%% Create the Phi and Gamma matrices
Phi = condensed_initial_gen(A, N);
Gamma = condensed_prediction_gen(A, B, N);


%% Create the final matrix
J = Gamma'*Qbar*Phi + Sbar'*Phi;

end

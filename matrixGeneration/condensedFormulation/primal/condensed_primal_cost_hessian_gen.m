function [ H ] = condensed_primal_cost_hessian_gen( N, A, B, Q, R, varargin )
%CONDENSED_PRIMAL_COST_HESSIAN_GEN Generate the Hessian matrix H
%
% Create the condensed Hessian matrix for the cost function of the condensed
% linear time-invariant MPC problem.
%
% If no P matrix is supplied, it defaults to Q.
% If no S matrix is supplied, it defaults to 0.
%
% Usage:
%   [ H ] = CONDENSED_PRIMAL_COST_HESSIAN_GEN( N, A, B, Q, R )
%   [ H ] = CONDENSED_PRIMAL_COST_HESSIAN_GEN( N, A, B, Q, R, P )
%   [ H ] = CONDENSED_PRIMAL_COST_HESSIAN_GEN( N, A, B, Q, R, P, S )
%
% Inputs:
%   N - The horizon length
%   A - The state transition matrix
%   B - The input mapping matrix
%   Q - The state weighting matrix
%   R - The input weighting matrix
%   P - The final state weighting matrix
%   S - The state-input cross term weight matrix
%
% Outputs:
%   H - The Hessian matrix
%
% See also CONDENSED_PRIMAL_COST_LINEAR_GEN
%
% Created by: Ian McInerney
% Created on: May 21, 2018
% Version: 1.2
% Last Modified: August 17, 2018
%
% Revision History
%   1.0 - Initial release  
%   1.1 - Added final state cost as separate input
%   1.2 - Added S term and rearranged inputs

[n, m] = size(B);

%% Parse the input arguments
p = inputParser;
addOptional(p, 'P', Q);
addOptional(p, 'S', zeros(n,m));
parse(p,varargin{:});

% Extract the matrices
P = p.Results.P;
S = p.Results.S;

% See if P was provided by the user
if (isempty(P))
    P = Q;
end


%% Create the prediction matrix gamma
gamma = condensed_prediction_gen(A, B, N);


%% Create the diagonal matrices for the Q, R and S matrices
Qbar = kron(eye(N-1), Q);
Qbar = blkdiag(Qbar, P);

Sbar = kron(eye(N-1), S);
Sbar = blkdiag(Sbar, zeros(n,m));

Rbar = kron(eye(N), R);


%% Put them all together
H = (gamma'*Qbar*gamma + gamma'*Sbar + Sbar'*gamma + Rbar);

end


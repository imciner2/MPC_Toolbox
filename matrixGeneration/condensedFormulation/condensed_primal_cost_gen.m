function [ H, J ] = condensed_primal_cost_gen( N, A, B, Q, R, varargin )
%CONDENSED_PRIMAL_COST_GEN Generate the matrices for the condensed primal cost
%
% Create the Hessian matrix and linear-term matrix for the cost function 
% of the condensed primal form for linear time-invariant MPC problem.
%
% Intial optimization problem:
%   min  x_n'Px_n + sum( x_k'Qx_k + 2x_k'Su_k + u_k'Ru_k )
%   s.t. Dx_k <= Cd
%        Eu_k <= Ce
%
% Resulting optimization problem:
%   min  u'Hu + x0'Ju
%   s.t. Gu <= Fx0 + g
%
% If no P matrix is supplied, it defaults to Q.
% If no S matrix is supplied, it defaults to 0.
%
%
% Usage:
%   [ H, J ] = CONDENSED_PRIMAL_COST_GEN( N, A, B, Q, R )
%   [ H, J ] = CONDENSED_PRIMAL_COST_GEN( N, A, B, Q, R, P )
%   [ H, J ] = CONDENSED_PRIMAL_COST_GEN( N, A, B, Q, R, P, S )
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
%   J - The linear-term matrix
%
% see also: CONDENSED_PRIMAL_CONSTRAINT_GEN
%
% Created by: Ian McInerney
% Created on: May 21, 2018
% Version: 1.3
% Last Modified: November 16, 2018
%
% Revision History
%   1.0 - Initial release  
%   1.1 - Added final state cost as separate input
%   1.2 - Added S term and rearranged inputs
%   1.3 - Added linear term

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


%% Create the prediction matrix and initial condition matrix
Phi   = condensed_initial_gen(A, N);
Gamma = condensed_prediction_gen(A, B, N);


%% Create the diagonal matrices for the Q, R and S matrices
Qbar = kron(eye(N-1), Q);
Qbar = blkdiag(Qbar, P);

Sbar = kron(eye(N-1), S);
Sbar = blkdiag(Sbar, zeros(n,m));

Rbar = kron(eye(N), R);


%% Form the Hessian
H = (Gamma'*Qbar*Gamma + Gamma'*Sbar + Sbar'*Gamma + Rbar);


%% Form the linear term
J = Gamma'*Qbar*Phi + Sbar'*Phi;

end

function [ G, F, g ] = gen_con_ineq( N, A, B, E, Ce, varargin )
%GEN_CON_INEQ Generate the inequality constraint matrices for the condensed formulation
%
% Create the condensed matrices for the inequality constraints of the
% condensed linear time-invariant MPC problem.
%
% Intial optimization problem:
%   min  0.5*x_n'*P*x_n + 0.5*sum( x_k'*Q*x_k + 2*x_k'*S*u_k + u_k'*R*u_k )
%   s.t. D*x_k <= Cd
%        E*u_k <= Ce
%
% Resulting optimization problem:
%   min  0.5*u'*H*u + x0'*J*u
%   s.t. G*u <= F*x0 + g
%
%
% Usage:
%   [ G, F, g ] = GEN_CON_INEQ( N, A, B, E, Ce )
%   [ G, F, g ] = GEN_CON_INEQ( N, A, B, E, Ce, D, Cd )
%
% Inputs:
%   N  - The horizon length
%   A  - The state transition matrix
%   B  - The input mapping matrix
%   E  - The stage constraints for the inputs
%   Ce - The right-hand side of the input constraints
%   D  - The stage constraints for the states
%   Cd - The right-hand side of the state constraints
%
% Outputs:
%   G - The coefficient matrix
%   F - The matrix multiplying the initial state
%   g - The constant vector
%
% See also CONDENSED.PRIMAL.GEN_COST
%
% Created by: Ian McInerney
% Created on: August 17, 2018
% Version: 1.1
% Last Modified: November 17, 2018
%
% Revision History
%   1.0 - Initial release
%   1.1 - Added remaining constraint matrices


%% Parse the input arguments
p = inputParser;
addOptional(p,  'D', []);
addOptional(p, 'Cd', []);
parse(p,varargin{:});

% Extract the matrices
D  = p.Results.D;
Cd = p.Results.Cd;


%% Find the number of constraints and the system size
[nE, ~] = size(E);
[nD, ~] = size(D);
[n, m]  = size(B);


%% Get the condensed system matrices
Phi   = condensed.gen_initial(A, N);
Gamma = condensed.gen_prediction(A, B, N);


%% Create the component matrices
Dbar = kron( eye(N), [D; zeros(nE, n)]);
Ebar = kron( eye(N), [zeros(nD,m); E]);


%% Compute the coefficient matrix
G = Dbar*Gamma + Ebar;


%% Compute the initial state matrix
F = -Dbar*Phi;


%% Compute the constant vector
g = kron( ones(N,1), [Cd; Ce] );

end

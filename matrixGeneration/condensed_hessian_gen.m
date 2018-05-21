function [ H ] = condensed_hessian_gen( A, B, Q, R, N )
%CONDENSED_HESSIAN_GEN Generate the condensed Hessian matrix
%
% Create the condensed Hessian matrix for the linear time-invariant MPC
% problem.
%
%
% Usage:
%   [ H ] = CONDENSED_HESSIAN_GEN( A, B, Q, R, N )
%
% Inputs:
%   A - The state transition matrix
%   B - The input mapping matrix
%   Q - The state weighting matrix
%   R - The input weighting matrix
%   N - The horizon length
%
% Outputs:
%   H - The Hessian matrix
%
%
% Created by: Ian McInerney
% Created on: May 21, 2018
% Version: 1.0
% Last Modified: May 21, 2018
%
% Revision History
%   1.0 - Initial release  


%% Create the prediction matrix gamma
gamma = condensed_prediction_gen(A, B, N);


%% Create the diagonal matrices for the Q and R matrices
Qbar = kron(eye(N+1), Q);
Rbar = kron(eye(N), R);


%% Put them all together
H = (gamma'*Qbar*gamma + Rbar);

end


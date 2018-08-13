function [ G ] = condensed_linear_gen( A, B, Q, P, N )
%CONDENSED_LINEAR_GEN Generate the condensed matrix for the linear term
%
% Create the condensed matrix for the linear term in the linear
% time-invariant MPC problem.
%
%
% Usage:
%   [ G ] = CONDENSED_LINEAR_GEN( A, B, Q, N )
%
% Inputs:
%   A - The state transition matrix
%   B - The input mapping matrix
%   Q - The state weighting matrix
%   P - The final state weighting matrix
%   N - The horizon length
%
% Outputs:
%   G - The linear term matrix
%
%
% Created by: Ian McInerney
% Created on: May 22, 2018
% Version: 1.0
% Last Modified: May 22, 2018
%
% Revision History
%   1.0 - Initial release  


%% Create the diagonal matrix for the Q matrix
Qbar = kron(eye(N-1), Q);
Qbar = blkdiag(Qbar, P);


%% Create the Phi and Gamma matrices
Phi = condensed_initial_gen(A, N);
Gamma = condensed_prediction_gen(A, B, N);


%% Create the final matrix
G = Gamma'*Qbar*Phi;

end
function [ H ] = gen_prediction( A, B, N )
%GEN_PREDICTION Generate the full condensed prediction matrix
%
% Generate the prediction matrix for a MPC problem that removes the state
% variables from the optimization problem and only contains the possible
% inputs.
%
% The variable vector has the form of:
%   [ u_0, u_1, u_2, ..., u_(N-1) ]
%
% This prediction matrix formulation is described in Section 8.2 of
%   F. Borrelli, A. Bemporad, and M. Morari, Predictive Control for Linear
%   And Hybrid Systems. Cambridge, UK: Cambridge University Press, 2017.
%
%
% Usage:
%   [ H ] = GEN_PREDICTION( A, B, N );
%
% Inputs:
%   A - The discrete-time system state transition matrix
%   B - The discrete-time system input matrix
%   N - The horizon length
%
% Outputs:
%   H - The prediction matrix (in sparse form)
%
%
% Created by: Ian McInerney
% Created on: January 18, 2018
% Version: 1.1
% Last Modified: June 7, 2018
%
% Revision History
%   1.0 - Initial release
%   1.1 - Sped up the computation process


%% Figure out the size of the system
[numStates, numInputs] = size(B);


%% Create the main diagonal of the matrix
H = speye(N);
H = kron(H, B);


%% Create the state matrices
S = speye(N);
S = kron(S, eye(numStates));

d = ones(N-1, 1);   % Create the appropriate length diagonal vector
D = diag(d, -1);    % Put 1s on the appropriate sub diagonal
D = kron(D, -A);    % Put the component onto the diagonal
S = S + D;


%% Find the final matrix
H = S\H;

end


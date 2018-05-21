function [ H ] = gen_pred_sparse_initial( A, B, N )
%GEN_PRED_SPARSE_INITIAL Generate prediction matrix with initial state as a
%variable
%
% Generate the prediction matrix for a MPC problem that assumes the initial
% state is included in the variable vector.
%
% The variable vector has the form of:
%   [ x_0, u_0, x_1, u_1, ..., u_(N-1), x_N ]
%
% This prediction matrix formulation is described inside of 
%   E. N. Hartley, J. L. Jerez, A. Suardi, J. M. Maciejowski,
%   E. C. Kerrigan, and G. A. Constantinides, “Predictive Control Using
%   an FPGA With Application to Aircraft Control,” IEEE Trans. Control
%   Syst. Technol., vol. 22, no. 3, pp. 1006–1017, 2014.
%
% Usage:
%   [ H ] = gen_pred_sparse_initial( A, B, N );
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
% Created on: January 17, 2018
% Version: 1.0
% Last Modified: January 17, 2018
%
% Revision History
%   1.0 - Initial release


%% Determine the matrix sizes
[numStates, ~] = size(A);
[~, numInputs] = size(B);

%% Create the matrix block to use
I = eye(numStates);
Z = zeros(numStates, numInputs);
comp1 = [-I, Z];
comp2 = [ A, B];


%% Create the main part of the Matrix
H = speye(N+1);
H = kron(H,comp1);

%% Remove the last columns
H = H(:, 1:(end-numInputs));

%% Add in the system dynamics
G = speye(N);
G = kron(G, comp2);

[~, c] = size(H);
Z1 = zeros( numStates, c);
Z2 = zeros(numStates*N, numStates);
G = [Z1;
      G, Z2];
H = H + G;



end


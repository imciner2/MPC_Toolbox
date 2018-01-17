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
%   Hartley, E. N., et. al. (2014). Predictive Control Using an FPGA With
%   Application to Aircraft Control. IEEE Transactions on Control Systems
%   Technology, 22(3), 1006–1017.
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
component = [-I, Z;
              A, B];


%% Create the main part of the Matrix
H = speye(N);
H = kron(H, component);


%% Add on the last columns of the matrix
Z = zeros(N*2*numStates - numStates, numStates);
comp = [Z;
       -I];
H = [ H, comp];

end


function [ H ] = condensed_prediction_gen( A, B, N )
%CONDENSED_PREDICTION_GEN Generate the full condensed prediction matrix
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
%   [ H ] = CONDENSED_PREDICTION_GEN( A, B, N );
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
% Version: 1.0
% Last Modified: January 18, 2018
%
% Revision History
%   1.0 - Initial release


%% Figure out the size of the system
[numStates, numInputs] = size(B);


%% Create the main diagonal of the matrix
H = speye(N);
H = kron(H, B);


%% Create the diagonal elements of the matrix
for i=1:1:(N-1)
    comp = A^(i)*B;     % Create the element
    
    d = ones(N-i, 1);   % Create the appropriate length diagonal vector
    D = diag(d, -i);    % Put 1s on the appropriate sub diagonal
    
    D = kron(D, comp);  % Put the component onto the diagonal
    
    H = H + D;          % Add the new matrix to H
end


%% Create the first row of all zeros
H = [sparse(numStates, numInputs*N);
     H];

end


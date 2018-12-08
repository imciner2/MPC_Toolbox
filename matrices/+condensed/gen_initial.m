function [ Phi ] = gen_initial( A, N )
%GEN_INITIAL Generate the initial state mapping matrix
%
% Create the condensed initial state mapping matrix for the linear
% time-invariant MPC problem.
%
%
% Usage:
%   [ Phi ] = GEN_INITIAL( A, N )
%
% Inputs:
%   A - The state transition matrix
%   N - The horizon length
%
% Outputs:
%   Phi - The mapping matrix
%
%
% Created by: Ian McInerney
% Created on: May 22, 2018
% Version: 1.0
% Last Modified: May 22, 2018
%
% Revision History
%   1.0 - Initial release  


%% Create some variables
inter = A;
Phi = [];


%% Iterate over the horizon creating the powers of A
for (i = 1:1:N)
    % Append the most recent power to the matrix
    Phi = [Phi;
           inter];

    inter = A*inter;
end

end
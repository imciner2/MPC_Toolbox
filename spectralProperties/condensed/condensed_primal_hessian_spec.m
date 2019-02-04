function [ eigen ] = condensed_primal_hessian_spec( A, B, C, Q, R, N )
%CONDENSED_PRIMAL_HESSIAN_SPEC Estimate the eigenvalues of the condensed primal Hessian
%
% Estimate the eigenvalues of the condensed hessian matrix for a linear
% time-invariant system as the horizon length goes to infinity. This
% function requires that the system is Schur stable (e.g. all eigenvalues
% of A are inside the unit circle). It utilizes a circulant approximation
% to the Hessian matrix, and samples it according to the desired horizon
% length.
%
%
% Usage:
%   [ eigen ] = CONDENSED_PRIMAL_HESSIAN_SPEC( A, B, C, Q, R, N )
%
% Inputs:
%   A - The state transition matrix
%   B - The input mapping matrix
%   C - The output mapping matrix
%   Q - The state weighting matrix
%   R - The input weighting matrix
%   N - The horizon length
%
% Outputs:
%   eigen - The estimated eigenvalues
%
%
% Created by: Ian McInerney
% Created on: May 21, 2018
% Version: 1.0
% Last Modified: May 21, 2018
%
% Revision History
%   1.0 - Initial release  


%% Create some helper variables
[n, m] = size(B);
I = eye(n);

% Figure out how many eigenvalues are needed per sample
numEigen = m;


%% Approximate the eigenvalues of the Hessian matrix
% Iterate over to approximate them
eigen = [];
for i = 0:1:(N-1)
    % Compute the complex number to use
    z = exp(1j*(-pi/2 + 2*pi*i/N));

    % Compute the matrix symbol
    M_c = B'*inv(z*I - A)'*C'*Q*C*inv(z*I - A)*B + R;
    
    % Compute the eigenvalues of the matrix symbol
    % abs is only here to prevent warnings, the eigenvalues should be
    % positive real anyway since the Hessian is real-symmetric
    e = abs( eigs(M_c, numEigen) );
    eigen = [eigen;
             e];
end
sort(eigen);

end


function [ k, varargout ] = condensed_hessian_cond( A, B, C, Q, R )
%CONDENSED_HESSIAN_COND Estimate the asymptotic condition number
%
% Estimate the asymptotic condition number of the condensed hessian matrix
% for a linear time-invariant system as the horizon length goes to
% infinity. This function requires that the system is Schur stable (e.g.
% all eigenvalues of A are inside the unit circle).
%
%
% Usage:
%   [ k ] = CONDENSED_HESSIAN_COND( A, B, C, Q, R )
%   [ k, maxE ] = CONDENSED_HESSIAN_COND( A, B, C, Q, R )
%   [ k, maxE, minE ] = CONDENSED_HESSIAN_COND( A, B, C, Q, R )
%
% Inputs:
%   A - The state transition matrix
%   B - The input mapping matrix
%   C - The output mapping matrix
%   Q - The state weighting matrix
%   R - The input weighting matrix
%
% Outputs:
%   k    - The condition number
%   maxE - The largest eigenvalue used in the computation
%   minE - The smallest eigenvalue used in the computation
%
%
% Created by: Ian McInerney
% Created on: May 21, 2018
% Version: 1.0
% Last Modified: May 21, 2018
%
% Revision History
%   1.0 - Initial release  


if ( (nargout-1) > 2 )
    error('Too many outputs requested');
    return
end

%% Create some helper variables
[n, m] = size(B);
I = eye(n);


%% Find the smallest and largest eigenvalues over a set range
num = 10000;
maxE = NaN;
minE = NaN;
for i = 0:1:(num-1)
    z = exp(1j*(-pi/2 + 2*pi*i/num));

    % Compute the matrix symbol
    M_c = B'*inv(z*I - A)'*Q*inv(z*I - A)*B + R;
    
    % Compute the eigenvalues of the matrix symbol
    % abs is only here to prevent warnings, the eigenvalues should be
    % positive real anyway
    me = abs( eigs(M_c, 1, 'LM') );
    maxE = max( [maxE, me]);
    
    me = abs( eigs(M_c, 1, 'SM') );
    minE = min( [minE, me]);
end


%% Compute the condition number
if (minE <= 1e-10)
    warning('Hessian is singular');
    k = inf;
else
    k = maxE./minE;
end


%% Output the eigenvalues if desired
switch( (nargout-1) )
    case 2
        varargout{2} = minE;
        varargout{1} = maxE;
    case 1
        varargout{1} = maxE;
end

end


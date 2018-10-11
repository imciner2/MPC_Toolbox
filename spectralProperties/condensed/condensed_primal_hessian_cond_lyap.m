function [ k, varargout ] = condensed_primal_hessian_cond_lyap( sys, Q, R )
%CONDENSED_PRIMAL_HESSIAN_COND_LYAP Estimate the asymptotic condition number
%
% Estimate the asymptotic condition number of the condensed hessian matrix
% for a linear time-invariant system as the horizon length goes to
% infinity. This function requires that the system is Schur stable (e.g.
% all eigenvalues of A are inside the unit circle).
%
% It assumes that the terminal cost matrix (P) is derived as the solution
% to the discrete-time Lyapunov equation.
%
% Usage:
%   [ k ] = CONDENSED_PRIMAL_HESSIAN_COND_LYAP( sys, Q, R )
%   [ k, maxE ] = CONDENSED_PRIMAL_HESSIAN_COND_LYAP( sys, Q, R )
%   [ k, maxE, minE ] = CONDENSED_PRIMAL_HESSIAN_COND_LYAP( sys, Q, R )
%
% Inputs:
%   sys - The discrete-time system being predicted
%   Q - The state weighting matrix
%   R - The input weighting matrix
%
% Outputs:
%   k    - The upper bound on the condition number
%   maxE - The upper bound on the smallest eigenvalue
%   minE - The lower bound on the largest eigenvalue
%
%
% See also: CONDENSED_PRIMAL_HESSIAN_COND_SAME
%
% Created by: Ian McInerney
% Created on: August 30, 2018
% Version: 1.0
% Last Modified: August 31, 2018
%
% Revision History
%   1.0 - Initial release


%% Extract the matrices
sys = ss(sys);
A = sys.A;
B = sys.B;
C = sys.C;
D = sys.D;


%% Create some helper variables
[n, m] = size(B);
I = eye(n);


%% Verify the output arguments requested
if ( (nargout-1) > 2 )
    error('Too many outputs requested');
    return
end


%% Get the eigenvalue estimates for the 
[~, maxE_q, minE_q] = condensed_primal_hessian_cond_same( sys, Q, R);


%% Compute the eigenvalues of the correction term Hp
% Compute the terminal weight
P = dlyap(A', Q);
P2 = sqrtm(P);

% Compute the observability gramian
W = dgram(sys.A, sys.B);

% Form the correction term and compute its eigenvalues
vv = P2*(W - sys.B*sys.B' )*P2;
e = eig( vv );


%% Compute the eigenvalue bounds
minE = minE_q;
maxE = maxE_q + max(e);


%% Compute the condition number bound
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


function [ k, varargout ] = condensed_primal_hessian_cond_same( sys, Q, R, varargin )
%CONDENSED_PRIMAL_HESSIAN_COND_SAME Estimate the asymptotic condition number
%
% Estimate the asymptotic condition number of the condensed hessian matrix
% for a linear time-invariant system as the horizon length goes to
% infinity. This function requires that the system is Schur stable (e.g.
% all eigenvalues of A are inside the unit circle).
%
% It assumes that the terminal cost weight matrix (P) is the same as the
% stage cost matrix (Q).
%
%
% Usage:
%   [ k ] = CONDENSED_PRIMAL_HESSIAN_COND_SAME( sys, Q, R )
%   [ k, maxE ] = CONDENSED_PRIMAL_HESSIAN_COND_SAME( sys, Q, R )
%   [ k, maxE, minE ] = CONDENSED_PRIMAL_HESSIAN_COND_SAME( sys, Q, R )
%   [ k, maxE, minE ] = CONDENSED_PRIMAL_HESSIAN_COND_SAME( sys, Q, R, S )
%
% Inputs:
%   sys - The discrete-time system being predicted
%   Q   - The state weighting matrix
%   R   - The input weighting matrix
%   S   - The state-input cross-term weighting matrix
%
% Outputs:
%   k    - The condition number
%   maxE - The largest eigenvalue
%   minE - The smallest eigenvalue
%
%
% See also: CONDENSED_PRIMAL_HESSIAN_COND_LYAP
%
% Created by: Ian McInerney
% Created on: May 21, 2018
% Version: 1.1
% Last Modified: August 31, 2018
%
% Revision History
%   1.0 - Initial release
%   1.1 - Added cross-term
%   1.2 - Modified to use system norms to speed up computation


%% Make sure it is a state-space system for easy access of the matrices
sys = ss(sys);


%% Create some helper variables
[n, m] = size( sys.B );
I = eye(n);


%% Make sure the D term is non-existent
if ( sum(sum( sys.D == 0 ) ) ~= n*m )
    warning('This function is not guaranteed to work on systems that contain a non-zero D matrix');
end


%% Make sure this is a discrete-time system
if ( sys.Ts == 0 )
    error('The dynamical system must be in discrete-time');
end


%% Parse the input arguments
p = inputParser;
addOptional(p, 'S', zeros(n,m));
parse(p,varargin{:});

% Extract the matrices
S = p.Results.S;


%% Verify the output arguments requested
if ( (nargout-1) > 2 )
    error('Too many outputs requested');
    return
end

%% Verify stability assumptions are present
if ( max( abs( eig(sys.A) ) ) >= 1)
    error('The given system is not Schur stable');
end


%% Create the Hessian matrix system
Hmsys = sys'*Q*sys + S'*sys + sys'*S + R;


%% Find the eigenvalue extremes
try
    % Utilize the H infinity norm to find the eigenvalues
    
    % Find the largest eigenvalue
    maxE = norm( Hmsys, 'inf' );
    
    % Find the smallest eigenvalue
    minE = norm( inv(Hmsys), 'inf' );
    if ( minE ~= 0 && minE ~= inf )
        minE = 1./minE;
    end

catch
    % Something broke when using the fast way, try the slow way now
    warning('Unable to use system norms, falling back to a global search');
    
    num = 10000;
    maxE = NaN;
    minE = NaN;
    for i = 0:1:(num-1)
        z = exp(1j*(-pi/2 + 2*pi*i/num));

        % Compute the inverse of the system matrix at this point
        tfm = evalfr( sys, z );

        % Compute the matrix symbol
        M_c = tfm'*Q*tfm + S'*tfm + tfm'*S + R;

        % Compute the eigenvalues of the matrix symbol
        % abs is only here to prevent warnings, the eigenvalues should be
        % positive real anyway
        me = abs( eigs(M_c, 1, 'LM') );
        maxE = max( [maxE, me]);

        me = abs( eigs(M_c, 1, 'SM') );
        minE = min( [minE, me]);
    end
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


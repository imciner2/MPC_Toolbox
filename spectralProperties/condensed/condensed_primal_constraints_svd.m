function [ minS, maxS ] = condensed_primal_constraints_svd( sys, E, varargin )
%CONDENSED_PRIMAL_CONSTRAINTS_SVD Estimate the extremal singular values
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
%   [ minS, maxS ] = CONDENSED_PRIMAL_CONSTRAINTS_SVD( sys, E)
%   [ minS, maxS ] = CONDENSED_PRIMAL_CONSTRAINTS_SVD( sys, E, D)
%
% Inputs:
%   sys - The discrete-time system being predicted
%   E   - The state weighting matrix
%   D   - The input weighting matrix
%
% Outputs:
%   minS - The smallest singular value
%   maxS - The largest singular value
%
%
% See also: CONDENSED_PRIMAL_HESSIAN_COND_LYAP
%
% Created by: Ian McInerney
% Created on: November 8, 2018
% Version: 1.0
% Last Modified: November 8, 2018
%
% Revision History
%   1.0 - Initial release


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
addOptional(p, 'D', []);
parse(p,varargin{:});

% Extract the matrices
D = p.Results.D;


%% Verify the output arguments requested
if ( (nargout-1) > 2 )
    error('Too many outputs requested');
    return
end

%% Verify stability assumptions are present
if ( max( abs( eig(sys.A) ) ) >= 1)
    error('The given system is not Schur stable');
end


%% Extract constraint information
[nic, ~] = size(E);
[nsc, ~] = size(D);


%% Create the matrix symbol
z = tf('z', sys.Ts);
sys1 = z*sys;

ieD = isempty(D);
ieE = isempty(E);

if ( ieD && ~ieE )
    % Only the inputs are present
    s = svd(E);
    
    maxS = max(s);
    minS = min(s);
    return;
    
elseif ( ~ieD && ieE )
    % Only the states are present
    Gsys = 1/z*D*sys.A*sys1 + D*sys.B;
    
elseif ( ~ieD && ~ieE )
    % Both matrices are present
    Dbar = [D*sys.A;
            zeros(nic, n)];
    Ebar = [D*sys.B;
            E];
    
    Gsys = Dbar*(1/z)*sys1 + Ebar;
    
else
    % Neither matrix is present
    warning('No constraint matrices provided');
    return
end


%% Find the extremal singular values
try
    % Utilize the H infinity norm to find the eigenvalues
    Gsysstar = Gsys'*Gsys;
    
    % Find the largest eigenvalue
    maxS = sqrt( norm( Gsysstar, 'inf', 1e-6 ) );
    
    % Find the smallest eigenvalue
    minS = sqrt( norm( inv(Gsysstar), 'inf', 1e-6 ) );
    if ( minS ~= 0 && minS ~= inf )
        minS = 1./minS;
    end

catch
    % Something broke when using the fast way, try the slow way now
    warning('Unable to use system norms, falling back to a global search');
    
    num = 30;
    maxS = NaN;
    minS = NaN;
    for i = 0:1:(num-1)
        z = exp(1j*(-pi/2 + 2*pi*i/num));

        % Compute the inverse of the system matrix at this point
        tfm = evalfr( sys1, z );

        % Compute the matrix symbol
        M_c = tfm'*Q*tfm + tfm'*S + S'*tfm + R;

        % Compute the eigenvalues of the matrix symbol
        % abs is only here to prevent warnings, the eigenvalues should be
        % positive real anyway
        me = abs( eigs(M_c, 1, 'LM') );
        maxS = max( [maxS, me]);

        me = abs( eigs(M_c, 1, 'SM') );
        minS = min( [minS, me]);
    end
end

end


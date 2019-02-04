function [ maxE, varargout ] = condensed_dual_hessian_eig( sys, Q, R, E, varargin )
%CONDENSED_DUAL_HESSIAN_EIG Compute the largest eigenvalue of the dual Hessian
%
% Usage:
%   [ maxE ] = CONDENSED_DUAL_HESSIAN_EIG( sys, Q, R, E );
%   [ maxE ] = CONDENSED_DUAL_HESSIAN_EIG( sys, Q, R, E, D );
%   [ maxE ] = CONDENSED_DUAL_HESSIAN_EIG( sys, Q, R, E, D, S );
%   [ maxE, minE ] = CONDENSED_DUAL_HESSIAN_EIG( sys, Q, R, E, D, opts );
%
% Inputs:
%   sys  - The physical system's model  
%   Q    - The Q matrix
%   R    - The R matrix
%   E    - The stage input constraints
%   D    - The stage state constraints
%   S    - The cross-term weight matrix
%   opts - Option pairs:
%          'Estimate' - 1 estimates the largest eigenvalue, 0 calculates it
%                       Defaults to 0
%          'Samples'  - The number of points to sample to find the largest
%                       eigenvalue. Defaults to 300
%
% Output:
%   maxE - The largest eigenvalue
%   minE - The smallest eigenvalue
%
%
% Created by: Ian McInerney
% Created on: November 18, 2018
% Version: 1.0
% Last Modified: November 18, 2018
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


%% Create the matrix symbol for the prediction matrix
z = tf('z', sys.Ts);
Pgam = z*sys;


%% Parse the input arguments
p = inputParser;
addOptional(p, 'D', []);
addOptional(p, 'S', []);
addParameter(p, 'Estimate', 0);
addParameter(p, 'Samples', 300);
parse(p,varargin{:});

% Extract the matrices
D = p.Results.D;
S = p.Results.S;
est = p.Results.Estimate;
nSamples = p.Results.Samples;


%% Determine the number of constraints
[nE, ~] = size(E);
[nD, ~] = size(D);


%% Create the matrix symbol for the constraints
Dbar = [D;
        zeros(nE, n)];

Ebar = [zeros(nD,m);
        E];

PG = Dbar*Pgam + Ebar;


%% If S is present, use the triangle-inequality bound
if ( ~isempty(S) )
    [~, ~, mH] = condensed_primal_hessian_cond_same(sys, Q, R, S);
    MG = sqrt( norm(PG'*PG, 'inf') );

    % If the calculation of the primal value failed, this one should also
    if (mH == inf)
        maxE = inf;
    else
        maxE = MG^2./mH;
    end
    
    return
end


%% Create the matrix symbol for the primal Hessian
PHp = Pgam'*Q*Pgam + R;


%% Create the matrix symbol for the dual Hessian's similar matrix
PHd1 = PG'*PG*inv(PHp);


%% Estimate the largest eigenvalue using the largest singular value
if ( est == 1 )
    maxE = norm( PHd1, 'inf' );
    return;
end


%% Find the largest eigenvalue
minE = NaN;
maxE = NaN;
for i = 0:1:(nSamples-1)
    z = exp(1j*(-pi/2 + 2*pi*i/nSamples));

    % Compute the matrix symbol at this point
    M_c = evalfr( PHd1, z );

    % Compute the eigenvalues of the matrix symbol
    % abs is only here to prevent warnings, the eigenvalues should be
    % positive real anyway
    me = abs( eigs(M_c, 1, 'LM') );
    maxE = max( [maxE, me]);
    
    mi = abs( eigs(M_c, 1, 'SM') );
    minE = min( [minE, mi]);
end

varargout{1} = minE;

end


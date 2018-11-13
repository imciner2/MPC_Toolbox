function [ maxE, varargout ] = condensed_dual_hessian_eig( sys, Q, R, E, varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

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
addParameter(p, 'Estimate', 0);
addParameter(p, 'Samples', 300);
parse(p,varargin{:});

% Extract the matrices
D = p.Results.D;
est = p.Results.Estimate;
nSamples = p.Results.Samples;


%% Determine the number of constraints
[nE, ~] = size(E);
[nD, ~] = size(D);


%% Create the matrix symbol for the primal Hessian
z = tf('z', sys.Ts);
Pgam = z*sys;
PHp = Pgam'*Q*Pgam + R;


%% Create the matrix symbol for the dual Hessian's similar matrix
Dbar = [D;
        zeros(nE, n)];

Ebar = [zeros(nD,m);
        E];

PG = Dbar*sys + Ebar;

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


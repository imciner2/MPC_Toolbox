function [ e ] = condensed_dual_hessian_spec_estimate( sys, N, Q, R, E, varargin )
%CONDENSED_DUAL_HESSIAN_SPEC_ESTIMATE Estimate the spectrum of the condensed dual Hessian
%
% This function will estimate the spectrum of the finite-horizon dual Hessian
% matrix using the Toeplitz theory (without actually forming the matrix).
%
%
% Usage:
% 	[ e ] = CONDENSED_DUAL_HESSIAN_SPEC_ESTIMATE( sys, N, Q, R, E );
% 	[ e ] = CONDENSED_DUAL_HESSIAN_SPEC_ESTIMATE( sys, N, Q, R, E, D );
%
% Inputs:
%   sys  - The physical system's model  
%   Q    - The Q matrix
%   R    - The R matrix
%   E    - The stage input constraints
%   D    - The stage state constraints
%
% Output:
%   e - The eigenvalues in sorted order
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
parse(p,varargin{:});

% Extract the matrices
D = p.Results.D;


%% Determine the number of constraints
[nE, ~] = size(E);
[nD, ~] = size(D);


%% Create the matrix symbol for the primal Hessian
PHp = Pgam'*Q*Pgam + R;


%% Create the matrix symbol for the dual Hessian's similar matrix
Dbar = [D;
        zeros(nE, n)];

Ebar = [zeros(nD,m);
        E];

PG = Dbar*Pgam + Ebar;

PHd1 = PG'*PG*inv(PHp);


%% Find the eigenvalues at the sampled points
e = [];
for i = 0:1:(N-1)
    z = exp(1j*(-pi/2 + 2*pi*i/N));

    % Compute the matrix symbol at this point
    M_c = evalfr( PHd1, z );

    % Compute the eigenvalues of the matrix symbol
    % abs is only here to prevent warnings, the eigenvalues should be
    % positive real anyway
    ei = abs( eig( M_c ) );
    e = [e;
         ei];
end


%% Sort the eigenvalues into numerical order
e = sort( e );

end


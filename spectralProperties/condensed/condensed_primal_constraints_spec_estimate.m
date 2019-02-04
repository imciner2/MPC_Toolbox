function [ s ] = condensed_primal_constraints_spec_estimate( sys, N, E, varargin )
%CONDENSED_PRIMAL_CONSTRAINTS_SPEC_ESTIMATE Estimate the spectrum of the condensed constraint matrix
%
% This function will estimate the singular value spectrum of the finite-horizon
% primal constraint matrix using the Toeplitz theory (without actually forming
% the matrix).
%
%
% Usage:
% 	[ s ] = CONDENSED_PRIMAL_CONSTRAINTS_SPEC_ESTIMATE( sys, N, E );
% 	[ s ] = CONDENSED_PRIMAL_CONSTRAINTS_SPEC_ESTIMATE( sys, N, E, D );
%
% Inputs:
%   sys  - The physical system's model  
%   E    - The stage input constraints
%   D    - The stage state constraints
%
% Output:
%   s - The singular values in sorted order
%
%
% Created by: Ian McInerney
% Created on: November 13, 2018
% Version: 1.0
% Last Modified: November 13, 2018
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
z = tf('z', sys.Ts);


%% Parse the input arguments
p = inputParser;
addOptional(p, 'D', []);
parse(p,varargin{:});

% Extract the matrices
D = p.Results.D;

% Extract the number of constraints
[nic, ~] = size(E);
[nsc, ~] = size(D);


%% Create the matrix symbol for the prediction matrix
z = tf('z', sys.Ts);
Pgam = z*sys;


%% Create the matrix symbol for the constraints
Dbar = [D;
        zeros(nE, n)];

Ebar = [zeros(nD,m);
        E];

PG = Dbar*Pgam + Ebar;


%% Find the singular values at the sampled points
s = [];
for i = 0:1:(N-1)
    z = exp(1j*(-pi/2 + 2*pi*i/N));

    % Compute the matrix symbol at this point
    if ( isempty(D) )
        M_c = PG;
    else
        M_c = evalfr( PG, z );
    end

    % Compute the singular values of the matrix symbol
    si = svd( M_c );
    s = [s;
         si];
end


% Sort the singular values into numerical order
s = sort( s );

end


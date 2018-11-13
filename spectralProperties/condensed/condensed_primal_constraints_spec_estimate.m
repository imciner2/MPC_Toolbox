function [ e ] = condensed_primal_constraints_spec_estimate( sys, N, E, varargin )
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


%% Create the Toeplitz symbol's system for the dual Hessian
sys1 = z*sys;

if ( isempty(D) )
    Gsys = E;
else
    Dsys = D*sys.A*(1/z)*sys1 + D*sys.B;
    
    if ( isempty(E) )
        Gsys = Dsys;
    else
        Gsys = [Dsys;
                 E];
    end
end


%% Find the largest eigenvalue
e = [];
for i = 0:1:(N-1)
    z = exp(1j*(-pi/2 + 2*pi*i/N));

    % Compute the matrix symbol at this point
    if ( isempty(D) )
        M_c = Gsys;
    else
        M_c = evalfr( Gsys, z );
    end

    % Compute the singular values of the matrix symbol
    ei = svd( M_c );
    e = [e;
         ei];
end

e = sort( e );

end


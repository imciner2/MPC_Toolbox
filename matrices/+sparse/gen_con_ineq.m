function [ Aineq, bineq ] = gen_con_ineq( N, n, m, E, Ce, varargin )
%GEN_CON_INEQ Generate the equality constraints for the linear dynamics
%
% Generate the prediction matrix for a MPC problem that assumes the initial
% state is included in the variable vector.
%
% The variable vector has the form of:
%   [ x_0, u_0, x_1, u_1, ..., u_(N-1), x_N ]
%
% This prediction matrix formulation is described inside of 
%   E. N. Hartley, J. L. Jerez, A. Suardi, J. M. Maciejowski,
%   E. C. Kerrigan, and G. A. Constantinides, “Predictive Control Using
%   an FPGA With Application to Aircraft Control,” IEEE Trans. Control
%   Syst. Technol., vol. 22, no. 3, pp. 1006–1017, 2014.
%
% Usage:
%   [ Aeq ] = GEN_CON_INEQ( N, n, m, E, Ce, n );
%   [ Aeq ] = GEN_CON_INEQ( N, n, m, E, Ce, D, Cd );
%
% Inputs:
%   N - The horizon length
%   n - The number of states
%   m - The number of inputs
%   D - The state stage constraints
%   E - The input stage constraints
%   Cd
%
% Outputs:
%   Aeq - The equality constraints for the linear dynamics
%
%
% Created by: Ian McInerney
% Created on: January 17, 2018
% Version: 1.0
% Last Modified: January 17, 2018
%
% Revision History
%   1.0 - Initial release


%% Parse the input arguments
p = inputParser;
addOptional(p,  'D', []);
addOptional(p, 'Cd', []);
parse(p,varargin{:});

% Extract the matrices
D  = p.Results.D;
Cd = p.Results.Cd;


%% Determine the variable sizes
[j, ~] = size(D);
[l, ~] = size(E);


%% Create the zero matrices
Zln = zeros(l, n);
Zjm = zeros(j, m);


%% Create the matrix for n=0 and n=N
if ( isempty(E) )
    A1 = [];
    b1 = [];
else
    A1 = [ Zln, E ];
    b1 = Ce;
end

if ( isempty(D) )
    Af = [];
    bf = [];
else
    Af = D;
    bf = Cd;
end


%% Create the matrix for the other stages
Ae = [   A1;
       Af, Zjm ];
be = [ b1;
       bf ];


%% Create the actual matrix
Aineq = speye(N-1);
Aineq = kron( Aineq, Ae );
Aineq = blkdiag( A1, Aineq, Af );

if ( isempty(Af) )
    [s, ~] = size(Aineq);
    Aineq = [Aineq, zeros( s, n)];
end


%% Create the RHS vectors
bineq = ones(N-1, 1);
bineq = kron( bineq, be );
bineq = [b1;
         bineq;
         bf];

end


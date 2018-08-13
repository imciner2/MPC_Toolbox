function [ Delta ] = fgm_coldStart_optim_delta( L, N, G, g, H, M )
%FGM_COLDSTART_OPTIM_DELTA Compute the optimal value for delta assuming cold-starting
%
% This function will compute the optimal value for Delta, which is used to
% create the lower iteration bound for the Fast Gradient Method.
% This requires the YALMIP toolbox to be installed, since it is used to
% interface to the optimization solvers
%
% This bound is derived in: 
%   S. Richter, C. N. Jones, and M. Morari, “Computational Complexity
%   Certification for Real-Time MPC With Input Constraints Based on the
%   Fast Gradient Method,” IEEE Transactions on Automatic Control,
%   vol. 57, no. 6, pp. 1391–1403, 2012.
%
%
% Usage:
%   [ delta ] = FGM_COLDSTART_OPTIM_DELTA( L, N, G, g);
%
% Inputs:
%   L    - The largest eigenvalue of the Hessian matrix
%   N    - The horizon length
%   G    - The input constraint matrix
%   g    - The input constraint vector
%
% Output:
%   delta - The delta parameter
%
%
% Created by: Ian McInerney
% Created on: August 13, 2018
% Version: 1.0
% Last Modified: August 13, 2018
%
% Revision History
%   1.0 - Initial release

%% Make sure thayt YALMIP is installed
if ( exist('yalmiptest', 'file') ~= 2 )
    error('YALMIP not found. Please install YALMIP and make sure it is on the path to run this function.');
end

%% Extract some variable sizing information
% r is the number of constraints in the constraint set
% m is the number of inputs in the constraint set
[r, m] = size(G);
[~, n] = size(M);
    
%% Create the larger matrices for G and g
Gl = kron( eye(N), G);
gl = kron( ones(N,1), g);


%% Create the variables
u = sdpvar(N*m,1);      % The optimal U variables
x = sdpvar(n, 1);       % State variables
lam = sdpvar(N*r, 1);   % The dual variables


%% The objective function
Obj = 0.5*u'*((L.*eye(N*m)) - H)*u;


%% Create the constraint functions
% Initial constraints
Con = [H*u + M*x + Gl'*lam == 0;
       Gl*u <= gl;
       0 <= lam];

% Create the complimentarity constraints
temp = (Gl*u - gl);
for i=1:1:(N*r)
    Con = [Con;
           lam(i)*temp(i) == 0];
end


%% Solve the optimization problem
ops = sdpsettings('verbose',0);
opt = optimize(Con, Obj, ops);

if (opt.problem ~= 0)
    yalErr = yalmiperror(opt.problem);
    error(['YALMIP error: ', yalErr]);
end


%% Extract the optimal cost, which is Delta
Delta = value(Obj);

end


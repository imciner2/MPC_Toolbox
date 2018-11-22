function [ y, i ] = dgp_quad( Hp, q, G, b, Ld, eps_g, varargin )
%DGP_QUAD Solve the quadratic program using the Dual Gradient Projection method
%
% This function uses the Dual Gradient Projection method (DGP) to solve a
% constrained quadratic program. The algorithm is based on the one
% presented in:
%   P. Patrinos, A. Guiggiani, and A. Bemporad, “A dual gradient-projection
%   algorithm for model predictive control in fixed-point arithmetic,”
%   Automatica, vol. 55, pp. 226–235, 2015.
%
% The problem takes in the primal matrices for the problem and computes the
% dual matrices itself as needed. Specifically, the primal problem is:
%   min  0.5*u'*Hp*u + q'*u
%   s.t. Gu <= b
% where for linear time-invariant MPC, q = J'*x0 and b = F*x0 + g.
% 
% This function returns the optimal dual solution. To recover the optimal
% primal solution, use the function CONDENSED_DUAL_PRIMAL_RECOVERY.
%
%
% Usage:
%   [ y ] = DGP_QUAD( Hp, q, G, b, Ld, eps_g )
%   [ y ] = DGP_QUAD( Hp, q, G, b, Ld, eps_g, iMax )
%   [ y, i ] = DGP_QUAD( Hp, q, G, b, Ld, eps_g )
%
% Inputs:
%   Hp    - The primal Hessian matrix
%   q     - The linear term of primal cost function
%   G     - The left-hand side of the primal inequality constraints
%   b     - The right-hand side of the primal inequality constraints
%   Ld    - The maximim eigenvalue of dual Hessian matrix
%   eps_g - The largest constraint satisfaction error allowed
%   iMax  - The Maximum number of iterations to do
%
% Outputs:
%   y - The dual variables
%   i - The number of iterations completed
%
%
% see also: CONDENSED_DUAL_PRIMAL_RECOVERY
%
% Created by: Ian McInerney
% Created on: November 18, 2018
% Version: 1.0
% Last Modified: November 18, 2018
%
% Revision History
%   1.0 - Initial release  


%% Find the size and create an identity matrix
[n, ~] = size(Hp);
I = eye(n);

[np, ~] = size(Hp);      % The number of primal variables
[nd, ~] = size(G);      % The number of dual variables


%% Parse the input arguments
p = inputParser;
addOptional(p, 'iMax', 1e6);
parse(p,varargin{:});

% Extract the inputs
iMax = p.Results.iMax;


%% Create some matrices for the computation
E = -(Hp\G');
e = -(Hp\q);


%% Create the initial condition
y  = zeros(nd, 1);
zs = zeros(np, 1);


%% Iterate until termination criteria is met
i = -1;     % Make -1 so that it immediately updates to 0
STOP = 0;
while ~STOP
    % Update the iteration count
    i = i+1;
    
    % Compute 
    z = E*y + e;
    
    % Compute the dual gradient
    gradY = G*z - b;
    
    % Update the dual variables in the new direction
    y = y + 1/Ld*gradY;
    y( y < 0 ) = 0;
    
    % Keep track of the primal iterates
    zs = zs + z;
    
    % Terminate after the iteration count maximum
    if ( i >= iMax )
        STOP = 1;
    end
    
    % Compute the stopping criteria
    zBar = 1/(i+1) * zs;
    er = G*zBar - b;
    er( er < 0 ) = 0;
    if ( norm(er, 'inf') < eps_g )
        STOP = 1;
    end
    
end


%% Return the number of iterations actually completed (not the iteration index)
i = i+1;

end

function [ k ] = fgm_iteration_bound( cond, L, eps, A, b )
%FGM_ITERATION_BOUND Compute the iteration bound for the FGM method applied
%to MPC
%
% Compute the lower iteration bound for the Fast Gradient Method when it is
% applied to a linear Model Predictive Control problem with input constraints
% of the form Au <= b.
%
% This lower bound is derived in 
%   S. Richter, C. N. Jones, and M. Morari, “Computational Complexity
%   Certification for Real-Time MPC With Input Constraints Based on the
%   Fast Gradient Method,” IEEE Transactions on Automatic Control,
%   vol. 57, no. 6, pp. 1391–1403, 2012.
%
%
% Usage:
%   [ k ] = FGM_ITERATION_BOUND( cond, L, eps, A, b );
%
% Inputs:
%   cond - The condition number of the Hessian matrix
%   L    - The largest eigenvalue of the Hessian matrix
%   eps  - The suboptimality level to solve to
%   A    - The input constraint matrix
%   b    - The input constraint vector
%
% Output:
%   k - The lower iteration bound
%
%
% Created by: Ian McInerney
% Created on: May 18, 2018
% Version: 1.0
% Last Modified: May 18, 2018
%
% Revision History
%   1.0 - Initial release  


%% Create the options for quadprog
opt = optimoptions('quadprog');
opt.Display = 'off';


%% Compute the Delta parameter
[n, m] = size(A);
I = eye(m);
Z = zeros(m, 1);

% Find the vector that has the largest 2-norm in the admissible input set
r = quadprog(I, Z, A, b, [], [], [], [], [], opt);
r = norm(r);

% Compute the conservative value of delta
delta = L/2 * r^2;


%% Compute two possible bounds
a1 = ceil( ( log(eps) - log(delta) ) / ( log(1 - sqrt(1/cond) ) ) );
a2 = ceil( 2*sqrt(delta/eps) - 2);


%% Choose the right bound
k = min( [a1, a2] );


%% Floor the iteration count to 0
k = max( [0, k] );

end


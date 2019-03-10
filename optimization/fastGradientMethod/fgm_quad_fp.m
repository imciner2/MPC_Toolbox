function [ u ] = fgm_quad_fp( H, b, projOp, beta, Linv, x0, maxIter )
%FGM_QUAD_FP Solve the quadratic program using the Fast Gradient Method
%
% This function uses the Fast Gradient Method (FGM) to solve a constrained 
% quadratic program. The constraints are handled through a projection
% operation, which is passed in as a function reference. The algorithm is
% based on the one presented in
%   J. L. Jerez, P. J. Goulart, S. Richter, G. A. Constantinides,
%   E. C. Kerrigan, and M. Morari, “Embedded Online Optimization for Model
%   Predictive Control at Megahertz Rates,” IEEE Transactions on Automatic
%   Control, vol. 59, no. 12, pp. 3238–3251, 2014.
%
% This algorithm solves an optimization problem of the form
%   min  0.5*u'*H*u + q'*u
%   s.t. G*u <= b
% where for linear time-invariant MPC, q = J'*x0 and b = F*x0 + g.
%
% The inequality constraint is handled through a projection operation onto
% the feasible set. The projection function is supplied by the user, and 
% takes in a single vector of variables, and then outputs a single vector
% of the variables after they are projected into the feasible set. For
% example, to implement box constraints -5 <= x <= 5, the projection
% function is:
%   projOp = @(x) min( 5, max(-5, x) );
%
% This algorithm is adapted for use with variables of the MATLAB embedded
% fixed-point type. This means it uses a constant step-size scheme and a
% termination criteria based on a maximum number of iterations.
%
%
% Usage:
%   [ u ] = FGM_QUAD_FP( H, b, projOp, beta, Linv, x0, maxIter )
%
% Inputs:
%   H       - The Hessian matrix
%   b       - The linear term of cost function
%   projOp  - A function to perform the projection operation
%   beta    - The step size to use
%   Linv    - The inverse of L
%   x0      - The intial point to start from
%   maxIter - The maximum number of iterations to perform
%
% Outputs:
%   u - The optimal control trajectory
%
%
% Created by: Ian McInerney
% Created on: March 10, 2019
% Version: 1.0
% Last Modified: March 10, 2019
%
% Revision History
%   1.0 - Initial release  


%% Find the size and create an identity matrix
[n, ~] = size(H);
I = eye(n);


%% Create the initial condition
z_n = x0;
y_n = x0;
y   = x0;
z   = x0;

%% Initialize some sizes to make the fixed-point types happy
gradY = b;
t = y;

%% Iterate until termination criteria is met
i = 0;
while i < maxIter
    i = i+1;
    
    % Copy the previous values into the current values
    y = y_n;
    z = z_n;
    
    % Compute the gradient
    gradY(:) = H*y + b;
    
    % Compute the new step
    t(:)   = y - Linv*gradY;
    
    % Project the new step onto the feasible space
    z_n(:) = projOp( t );
    
    % Use acceleration to get the next point
    y_n(:) = z_n + beta*(z_n - z);
end


%% Save the sequence for output
u = z_n;


end

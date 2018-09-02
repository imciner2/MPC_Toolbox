function [ u, i ] = fgm_quad( H, b, projOp, L, mu, varargin )
%FGM_QUAD Solve the quadratic program using the Fast Gradient Method
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
% The projection function takes in a single vector of variables, and then
% outputs a single vector of the variables after they are projected into
% the feasible space. For example, to implement box constraints -5 <= x <= 5,
% the projection function is
%   projOp = @(x) min( 5, max(-5, x) );
%
% This algorithm can utilize either a constant step size scheme where beta is
% held the same across all iterations (calculated using the formula in the 
% above paper), or it can utilize a variable stepsize.
%
% This algorithm can utilize one of 4 different stopping criteria:
%   'Iterations' - Terminate after the specified number of iterations
%   'Gradient'   - Terminate once the gradient falls below the given tolerance
%   'Conjugate'  - Terminate once the conjugate function falls below the tolerance
%   'Best'       - Combines Gradient and Conjugate to terminate when either
%                  is satisfied
% These stopping criteria are explained inside Chapter 6 of 
%   S. Richter, “Computational complexity certification of gradient
%   methods for real-time model predictive control,” PhD. Thesis,
%   ETH Zurich, 2012.
%
%
% Usage:
%   [ u ] = FGM_QUAD( H, b, projOp, L, mu, x0, 'stoppingCriteria', 'stoppingValue', 'stepsize' )
%
% Inputs:
%   H      - The Hessian matrix
%   b      - The linear term of cost function
%   projOp - A function to perform the projection operation
%   L      - The maximim eigenvalue of H
%   mu     - The minimum eigenvalue of H
%   x0     - The intial point to start from
%   stoppingCriteria - The method to use to stop the iterations
%   stoppingValue    - The value to stop the iterations at
%   stepsize         - Select either 'Variable' or 'Constant'
%
% Outputs:
%   u - The optimal control trajectory
%   i - The number of iterations completed
%
%
% Created by: Ian McInerney
% Created on: June 5, 2018
% Version: 1.1
% Last Modified: September 2, 2018
%
% Revision History
%   1.0 - Initial release  
%   1.1 - Added termination criteria


%% Find the size and create an identity matrix
[n, ~] = size(H);
I = eye(n);


%% Parse the input arguments
p = inputParser;
addRequired(p, 'x0');
addRequired(p, 'stoppingCriteria', @(x) isstring(x) || ischar(x));
addRequired(p, 'stoppingValue');
addRequired(p, 'stepsize', @(x) isstring(x) || ischar(x));
parse(p,varargin{:});

% Extract the inputs
x0 = p.Results.x0;
stoppingCriteria = p.Results.stoppingCriteria;
eps = p.Results.stoppingValue;
stepping = p.Results.stepsize;


%% Create the stopping criteria function
switch (stoppingCriteria)
    case 'Gradient'
        % Create a stopping criteria using the gradient method
        const = 0.5*(1/mu - 1/L);
        termFunc = @(i, x, y, x_n, y_n) (i == 0) || (const*norm(L*(y - x_n), 2)^2 > eps) ;
    case 'Conjugate'
        % Create a stopping criteria using the conjugacy method
        termFunc = @(i, x, y, x_n, y_n) x_n'*(H*x_n + b) + norm(H*x_n + b, 1) > eps;
    case 'Iterations'
        % Create an upper-bound using the iterations
        termFunc = @(i, x, y, x_n, y_n) i <= stoppingValue;
    case 'Best'
        % Combine the gradient method and conjugacy method
        const = 0.5*(1/mu - 1/L);
        termFunc = @(i, x, y, x_n, y_n) (i == 0) || ( (x_n'*(H*x_n + b) + norm(H*x_n + b, 1) > eps) && (const*norm(L*(y - x_n), 2)^2 > eps) );
    otherwise
        error('Unknown termination criteria');
end


%% Compute the step-size
switch(stepping)
    case 'Constant'
        beta = (sqrt(L) - sqrt(mu)) / (sqrt(L) + sqrt(mu));
    case 'Variable'
        alpha = sqrt(mu/L);
        alpha_n = sqrt(mu/L);
        stepping = 1;
end


%% Create the initial condition
x_n = x0;
y_n = x0;
y   = x0;
x   = x0;


%% Iterate until termination criteria is met
i = 0;
while termFunc(i, x, y, x_n, y_n)
    i = i+1;
    
    % Copy the previous values into the current values
    y = y_n;
    x = x_n;
    
    % Compute the gradient
    gradY = H*y + b;
    
    % Compute the new step
    t   = y - (1/L)*gradY;
    
    % Project the new step onto the feasible space
    x_n = projOp( t );
    
    % Compute a new value for alpha if desired
    if ( stepping == 1 )
        alpha = alpha_n;
        
        alphaFunc = @(x) (1-x)*alpha^2 + mu/L*x - x^2;
        alpha_n = fsolve( alphaFunc, alpha, optimoptions('fsolve','Display','off'));
        alpha_n = max( 0, min(alpha_n, 1));
        
        beta = (alpha*(1-alpha)) / (alpha^2 + alpha_n);
    end
    
    % Use acceleration to get the next point
    y_n = (1 + beta)*x_n - beta*x;
end


%% Save the sequence for output
u = x_n;


end

function [ u, varargout  ] = fgm_quad( H, b, projOp, L, mu, Imax )
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
%
% The projection function takes in a single vector of variables, and then
% outputs a single vector of the variables after they are projected into
% the feasible space. For example, to implement box constraints -5 <= x <= 5,
% the projection function is
%   projOp = @(x) min( 5, max(-5, x) );
%
%
% Usage:
%   [ u ] = FGM_QUAD( H, b, projOp, L, mu, Imax )
%   [ u, r ] = FGM_QUAD( H, b, projOp, L, mu, Imax )
%
% Inputs:
%   H      - The Hessian matrix
%   b      - The linear term of cost function
%   projOp - A function to perform the projection operation
%   L      - The maximim eigenvalue of H
%   mu     - The minimum eigenvalue of H
%   Imax   - The maximum number of iterations to run
%
% Outputs:
%   u - The optimal control trajectory
%   r - The 2-norm of the residuals from each iteration
%
%
% Created by: Ian McInerney
% Created on: June 5, 2018
% Version: 1.0
% Last Modified: June 5, 2018
%
% Revision History
%   1.0 - Initial release  

%% Save the residuals
if (nargout == 2)
    r = zeros(Imax, 1);
    saveRes = 1;
else
    saveRes = 0;
end


%% Find the size and create an identity matrix
[n, ~] = size(H);
I = eye(n);


%% Compute the step-size
beta = (sqrt(L) - sqrt(mu)) / (sqrt(L) + sqrt(mu));


%% Create the initial condition
z = ones(n,1);
y = z;


%% Iterate the desired number of times
for (i=1:1:(Imax+1))
    % Compute the new step
    t   = y - (1/L)*(H*y + b);
    
    if (saveRes)
        r(i) = norm(t, 2);
    end
    
    % Project the new step onto the feasible space
    z_n = projOp( t );
    
    % Use acceleration to get the next point
    y_n = (1 + beta)*z_n - beta*z;
    
    % Update the variables for the next iteration
    y = y_n;
    z = z_n;
end


%% Save the sequence for output
u = z;

if (saveRes)
    varargout{1} = r;
end

end